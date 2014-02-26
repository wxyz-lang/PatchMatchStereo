
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <vector>
#include <stack>
#include <list>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <omp.h>
#include "SLIC.h"
#include "Utilities.h"

#ifdef _DEBUG
#pragma comment(lib, "opencv_core248d.lib")
#pragma comment(lib, "opencv_highgui248d.lib")
#pragma comment(lib, "opencv_imgproc248d.lib")
#pragma comment(lib, "opencv_features2d248d.lib")
#pragma comment(lib, "opencv_calib3d248d.lib")
#pragma comment(lib, "opencv_video248d.lib")
#pragma comment(lib, "opencv_flann248d.lib")
#else
#pragma comment(lib, "opencv_core248.lib")
#pragma comment(lib, "opencv_highgui248.lib")
#pragma comment(lib, "opencv_imgproc248.lib")
#pragma comment(lib, "opencv_features2d248.lib")
#pragma comment(lib, "opencv_calib3d248.lib")
#pragma comment(lib, "opencv_video248.lib")
#pragma comment(lib, "opencv_flann248.lib")
#endif


#define USE_OPENMP
//#define LOAD_RESULT_FROM_LAST_RUN
#define DO_POST_PROCESSING
//#define USE_NELDERMEAD_OPT

// Static class member initialization 
std::stack<clock_t> Timer::time_stamps = std::stack<clock_t>();

// Gobal variables
int nrows, ncols;
const float BAD_PLANE_PENALTY = 120;  // defined as 2 times the max cost of dsi.
const float	gamma_proximity = 25;
int			g_improve_cnt = 0;

const int		patch_w = 35;
const int		patch_r = 17;
const int		maxiters = 3;
const float		alpha = 0.9;
const float		gamma = 10;
const float		tau_col = 10;
const float		tau_grad = 2;
const float		granularity = 0.25f;

const int folder_id = 3;
const std::string folders[] = { "tsukuba/", "venus/", "teddy/", "cones/", "Bowling2/", "Baby1/" };
const int scales[]	= { 16, 8, 4, 4, 3, 3 };
const int drange[]	= { 16, 20, 60, 60, 70, 70 };
const int scale				= scales[folder_id];
const int ndisps			= drange[folder_id];
const int dmax				= ndisps - 1;






inline bool InBound(float y, float x) { return 0 <= y && y < nrows && 0 <= x && x < ncols; }

VECBITMAP<float> ComputeColGradFeature(cv::Mat& img)
{
	int nrows = img.rows, ncols = img.cols;
	int sobel_scale = 1, sobel_delta = 0;
	cv::Mat gray, grad_x, grad_y;

	cv::cvtColor(img, gray, CV_BGR2GRAY);
	cv::Sobel(gray, grad_x, CV_32F, 1, 0, 3, sobel_scale, sobel_delta, cv::BORDER_DEFAULT);
	cv::Sobel(gray, grad_y, CV_32F, 0, 1, 3, sobel_scale, sobel_delta, cv::BORDER_DEFAULT);
	grad_x = grad_x / 8.f;
	grad_y = grad_y / 8.f;

	VECBITMAP<float> colgrad(nrows, ncols, 5);
	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			colgrad.get(y, x)[0] = img.at<cv::Vec3b>(y, x)[0];
			colgrad.get(y, x)[1] = img.at<cv::Vec3b>(y, x)[1];
			colgrad.get(y, x)[2] = img.at<cv::Vec3b>(y, x)[2];
			colgrad.get(y, x)[3] = grad_x.at<float>(y, x);
			colgrad.get(y, x)[4] = grad_y.at<float>(y, x);
		}
	}

	return colgrad;
}

VECBITMAP<float> ComputeAdGradientCostVolume(cv::Mat& imL, cv::Mat& imR, int ndisps, int sign, float granularity)
{
	int nrows = imL.rows, ncols = imL.cols;
	int nlevels = ndisps / granularity;

	VECBITMAP<float> colgradL = ComputeColGradFeature(imL);
	VECBITMAP<float> colgradR = ComputeColGradFeature(imR);
	VECBITMAP<float> dsiL(nrows, ncols, nlevels);

	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			for (int level = 0; level < nlevels; level++) {
				float d = level * granularity;
				float xm = x + sign * d;

				// FIXME: has implicitly assumed "ndisps <= ncols", it's not safe.
				if (xm < 0)			xm += ncols;
				if (xm > ncols - 1) xm -= ncols;

				float xmL, xmR, wL, wR;
				float cost_col, cost_grad;
				xmL = (int)(xm);
				xmR = (int)(xm + 0.99);
				wL = xmR - xm;
				wR = 1.f - wL;

				float *pL = colgradL.get(y, x);
				float *pRmL = colgradR.get(y, xmL);
				float *pRmR = colgradR.get(y, xmR);
				float pR[5];
				for (int i = 0; i < 5; i++) {
					pR[i] = wL * pRmL[i] + wR * pRmR[i];
				}

				cost_col = fabs(pL[0] - pR[0])
					+ fabs(pL[1] - pR[1])
					+ fabs(pL[2] - pR[2]);
				cost_col = std::min(tau_col, cost_col);
				cost_grad = fabs(pL[3] - pR[3])
					+ fabs(pL[4] - pR[4]);
				cost_grad = std::min(tau_grad, cost_grad);

				dsiL.get(y, x)[level] = (1 - alpha) * cost_col + alpha * cost_grad;
			}
		}
	}

	return dsiL;
}

unsigned hamdist(long long x, long long  y)
{
	int dist = 0;
	long long val = x ^ y;
	while (val) {
		++dist;
		val &= val - 1;
	}
	return dist;
}

VECBITMAP<long long> ComputeCensusImage(VECBITMAP<unsigned char> im)
{
	int vpad = 3, hpad = 4;
	VECBITMAP<long long> census(nrows, ncols);
	memset(census.data, 0, nrows * ncols * sizeof(long long));

	for (int yc = 0; yc < nrows; yc++) {
		for (int xc = 0; xc < ncols; xc++) {

			int uu = std::max(yc - vpad, 0);
			int dd = std::min(yc + vpad, nrows - 1);
			int ll = std::max(xc - hpad, 0);
			int rr = std::min(xc + hpad, ncols - 1);

			int idx = 0;
			long long feature = 0;
			unsigned char center = im[yc][xc];

			for (int y = uu; y <= dd; y++) {
				for (int x = ll; x <= rr; x++) {
					feature |= ((long long)(im[y][x] > center) << idx);
					idx++;
				}
			}

			census[yc][xc] = feature;
		}
	}

	return census;
}

VECBITMAP<float> ComputeCensusTensor(VECBITMAP<unsigned char>& imL, VECBITMAP<unsigned char>& imR, int sign)
{
	VECBITMAP<float> dsi(nrows, ncols, ndisps);
	VECBITMAP<long long> censusL = ComputeCensusImage(imL);
	VECBITMAP<long long> censusR = ComputeCensusImage(imR);

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			for (int d = 0; d < ndisps; d++) {
				int xm = (x + sign*d + ncols) % ncols;
				dsi.get(y, x)[d] = hamdist(censusL[y][x], censusR[y][xm]);
			}
		}
	}
	return dsi;
}

VECBITMAP<float> ComputeAdCensusCostVolume(cv::Mat& cvimL, cv::Mat& cvimR, int ndisps, int sign)
{
	const float ad_lambda = 30;
	const float census_labmda = 10;

	cv::Mat cvgrayL, cvgrayR;
	cv::cvtColor(cvimL, cvgrayL, CV_BGR2GRAY);
	cv::cvtColor(cvimR, cvgrayR, CV_BGR2GRAY);
	VECBITMAP<unsigned char> grayL(nrows, ncols, 1, cvgrayL.data);
	VECBITMAP<unsigned char> grayR(nrows, ncols, 1, cvgrayR.data);
	VECBITMAP<float> dsi_census = ComputeCensusTensor(grayL, grayR, sign);

	assert(cvimL.isContinuous());
	assert(cvimR.isContinuous());
	VECBITMAP<unsigned char> imL(nrows, ncols, 3, cvimL.data);
	VECBITMAP<unsigned char> imR(nrows, ncols, 3, cvimR.data);

	VECBITMAP<float> dsi(nrows, ncols, ndisps);
	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			for (int d = 0; d < ndisps; d++) {
				int xm = (x + sign*d + ncols) % ncols;
				unsigned char *pL = imL.get(y, x);
				unsigned char *pR = imR.get(y, xm);
				dsi.get(y, x)[d] = (fabs((float)pL[0] - pR[0])
					+ fabs((float)pL[1] - pR[1])
					+ fabs((float)pL[2] - pR[2])) / 3.f;
			}
		}
	}

	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			for (int d = 0; d < ndisps; d++) {
				dsi.get(y, x)[d] = 2 - exp(-dsi.get(y, x)[d] / ad_lambda) - exp(-dsi_census.get(y, x)[d] / census_labmda);
			}
		}
	}

	return dsi;
}

VECBITMAP<float> WinnerTakesAll(VECBITMAP<float>& dsi)
{
	int nrows = dsi.h, ncols = dsi.w, ndisps = dsi.n;
	VECBITMAP<float> disp(nrows, ncols);

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			int minidx = 0;
			float *cost = dsi.get(y, x);
			for (int k = 1; k < ndisps; k++) {
				if (cost[k] < cost[minidx]) {
					minidx = k;
				}
			}
			disp[y][x] = minidx;
		}
	}
	return disp;
}

void LocalSearch(cv::Mat& imL, cv::Mat& imR, int ndisps, cv::Mat& dispL, cv::Mat& dispR)
{
	VECBITMAP<float> dsiL = ComputeAdGradientCostVolume(imL, imR, ndisps, -1, 1.f);
	VECBITMAP<float> dsiR = ComputeAdGradientCostVolume(imR, imL, ndisps, +1, 1.f);

	VECBITMAP<float> dL = WinnerTakesAll(dsiL);
	VECBITMAP<float> dR = WinnerTakesAll(dsiR);

	int nrows(imL.rows), ncols(imL.cols);
	dispL.create(nrows, ncols, CV_32FC1);
	dispR.create(nrows, ncols, CV_32FC1);

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			dispL.at<float>(y, x) = dL[y][x];
			dispR.at<float>(y, x) = dR[y][x];
		}
	}
}

double ComputePlaneCost(int yc, int xc, Plane& coeff_try, VECBITMAP<float>& dsi, VECBITMAP<float>& w)
{
	double cost = 0;
	for (int y = yc - patch_r; y <= yc + patch_r; y++) {
		for (int x = xc - patch_r; x <= xc + patch_r; x++) {
			float d = (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
			int level = 0.5 + d / granularity;
			if (InBound(y, x)) {
				if (d < 0 || d > dmax) {	// must be a bad plane.
					cost += BAD_PLANE_PENALTY;
				}
				else {
					cost += w[y - yc + patch_r][x - xc + patch_r] * (dsi.get(y, x)[level]);
				}
			}
		}
	}
	return cost;
}

void RandomInit(VECBITMAP<Plane>& coeffs, VECBITMAP<float>& bestcosts, VECBITMAP<float>& dsi, VECBITMAP<float>& weights)
{
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			coeffs[y][x].RandomAssign(y, x, dmax);
		}
	}
	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			bestcosts[y][x] = ComputePlaneCost(y, x, coeffs[y][x], dsi, VECBITMAP<float>(patch_w, patch_w, 1, weights.line_n1(y*ncols + x)));
		}
	}
}

VECBITMAP<float> PrecomputeWeights(cv::Mat& img)
{
	assert(img.isContinuous());
	VECBITMAP<unsigned char> im(nrows, ncols, 3, img.data);

	const int nw = nrows * ncols * patch_w * patch_w;
	const int patchsize = patch_w * patch_w;

	VECBITMAP<float> ret(nrows * ncols, patch_w * patch_w);
	float *weights = ret.data;
	memset(weights, 0, nw * sizeof(float));

	float w_prox[patch_w * patch_w];
	for (int y = 0; y < patch_w; y++) {
		for (int x = 0; x < patch_w; x++) {
			float dist = sqrt((y - patch_r) * (y - patch_r) + (x - patch_r) * (x - patch_r));
			w_prox[y * patch_w + x] = exp(-dist / gamma_proximity);
		}
	}

	#pragma omp parallel for
	for (int yc = 0; yc < nrows; yc++) {
		for (int xc = 0; xc < ncols; xc++) {

			float *w = &weights[(yc*ncols + xc) * patchsize];
			unsigned char *rgb1 = im.get(yc, xc);
			unsigned char *rgb2 = NULL;

			int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
			int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

			for (int y = yb; y <= ye; y++) {
				for (int x = xb; x <= xe; x++) {
					rgb2 = im.get(y, x);
					float cost = std::abs((float)rgb1[0] - (float)rgb2[0])
						+ std::abs((float)rgb1[1] - (float)rgb2[1])
						+ std::abs((float)rgb1[2] - (float)rgb2[2]);
					w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] = cost;
				}
			}
		}
	}

	for (int i = 0; i < nw; i++) {
		weights[i] = exp(-weights[i] / gamma);
	}

	return ret;
}

void ImproveGuess(int y, int x, Plane& coeff_old, float& bestcost, Plane& coeff_try, VECBITMAP<float>& dsi, VECBITMAP<float>& w)
{
	float cost = ComputePlaneCost(y, x, coeff_try, dsi, w);
	if (cost < bestcost) {
		g_improve_cnt++;
		bestcost = cost;
		coeff_old = coeff_try;
	}
}

void PropagateAndRandomSearch(int y, int x,
	VECBITMAP<Plane>& coeffsL,		VECBITMAP<Plane>& coeffsR,
	VECBITMAP<float>& bestcostsL,	VECBITMAP<float>& bestcostsR,
	VECBITMAP<float>& dsiL,			VECBITMAP<float>& dsiR,
	VECBITMAP<float>& weightsL,		VECBITMAP<float>& weightsR,
	int iter, int sign)
{
	int xchange, ychange;
	if (iter % 2 == 0)  xchange = ychange = +1;
	else				xchange = ychange = -1;
	VECBITMAP<float> wL(patch_w, patch_w, 1, weightsL.line_n1(y*ncols + x));

#ifndef USE_NELDERMEAD_OPT
	// Spatial Propagation
	int qy = y - ychange, qx = x;
	if (InBound(qy, qx)) {
		Plane coeff_try = coeffsL[qy][qx];
		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
	}

	qy = y; qx = x - xchange;
	if (InBound(qy, qx)) {
		Plane coeff_try = coeffsL[qy][qx];
		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
	}

	// Random Search
	float radius_z = dmax / 2.0f;
	float radius_n = 1.0f;
	while (radius_z >= 0.1) {
		Plane coeff_try = coeffsL[y][x].RandomSearch(y, x, radius_z, radius_n, dmax);
		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
		radius_z /= 2.0f;
		radius_n /= 2.0f;
	}

	// View Propagation
	Plane coeff_try = coeffsL[y][x].ReparametrizeInOtherView(y, x, sign, qy, qx);
	if (0 <= qx && qx < ncols) {
		VECBITMAP<float> wR(patch_w, patch_w, 1, weightsR.line_n1(qy*ncols + qx));
		ImproveGuess(qy, qx, coeffsR[qy][qx], bestcostsR[qy][qx], coeff_try, dsiR, wR);
	}
#else
	const int dx[] = { -1, 0, +1, 0 };
	const int dy[] = { 0, -1, 0, +1 };
	float nm_opt_x[3 * 4];

	for (int dir = 0; dir < 4; dir++) {
		int qy = y + dy[dir];
		int qx = x + dx[dir];
		if (InBound(qy, qx)) {
			coeffsL[qy][qx].GetAbc(nm_opt_x + 3 * dir);
		}
		else {
			RandomAbc(nm_opt_x + 3 * dir, qy, qx);
		}
	}
	nm_opt_struct.yc = y;
	nm_opt_struct.xc = x;
	nm_opt_struct.w = &wL;
	nm_opt_struct.dsi = &dsiL;

	NelderMeadOptimize(nm_opt_x, nm_compute_plane_cost, 5);
	coeffsL[y][x].SetAbc(nm_opt_x);
#endif
}

void PlaneMapToDisparityMap(VECBITMAP<Plane>& coeffs, VECBITMAP<float>& disp)
{
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			disp[y][x] = coeffs[y][x].ToDisparity(y, x);
		}
	}
}

void CrossCheck(VECBITMAP<float>& dispL, VECBITMAP<float>& dispR, VECBITMAP<bool>& validL, VECBITMAP<bool>& validR)
{
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {

			int xR = std::max(0.f, std::min((float)ncols, x - dispL[y][x]));
			validL[y][x] = (std::abs(dispL[y][x] - dispR[y][xR]) <= 1);

			int xL = std::max(0.f, std::min((float)ncols, x + dispR[y][x]));
			validR[y][x] = (std::abs(dispR[y][x] - dispL[y][xL]) <= 1);
		}
	}
}

void FillHole(int y, int x, VECBITMAP<bool>& valid, VECBITMAP<Plane>& coeffs)
{
	// This function fills the invalid pixel (y,x) by finding its nearst (left and right) 
	// valid neighbors on the same scanline, and select the one with lower disparity.

	int xL = x - 1, xR = x + 1, bestx = x;
	while (!valid[y][xL] && 0 <= xL) {
		xL--;
	}
	while (!valid[y][xR] && xR < ncols) {
		xR++;
	}
	if (0 <= xL) {
		bestx = xL;
	}
	if (xR < ncols) {
		if (bestx == xL) {
			float dL = coeffs[y][xL].ToDisparity(y, x);
			float dR = coeffs[y][xR].ToDisparity(y, x);
			if (dR < dL) {
				bestx = xR;
			}
		}
		else {
			bestx = xR;
		}
	}
	coeffs[y][x] = coeffs[y][bestx];
}

void WeightedMedianFilter(int yc, int xc, VECBITMAP<float>& disp, VECBITMAP<float>& weights, VECBITMAP<bool>& valid, bool useInvalidPixels)
{
	std::vector<std::pair<float, float>> dw_pairs;

	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

	for (int y = yb; y <= ye; y++) {
		for (int x = xb; x <= xe; x++) {
			if (useInvalidPixels || valid[y][x]) {
				std::pair<float, float> dw(disp[y][x], weights[y - yc + patch_r][x - xc + patch_r]);
				dw_pairs.push_back(dw);
			}
		}
	}

	std::sort(dw_pairs.begin(), dw_pairs.end());

	float w = 0.f, wsum = 0.f;
	for (int i = 0; i < dw_pairs.size(); i++) {
		wsum += dw_pairs[i].second;
	}

	for (int i = 0; i < dw_pairs.size(); i++) {
		w += dw_pairs[i].second;
		if (w >= wsum / 2.f) {
			// Note that this line can always be reached.
			if (i > 0) {
				disp[yc][xc] = (dw_pairs[i - 1].first + dw_pairs[i].first) / 2.f;
			}
			else {
				disp[yc][xc] = dw_pairs[i].first;
			}
			break;
		}
	}
}

void PostProcess(
	VECBITMAP<float>& weightsL, VECBITMAP<float>& weightsR,
	VECBITMAP<Plane>& coeffsL,	VECBITMAP<Plane>& coeffsR,
	VECBITMAP<float>& dispL,	VECBITMAP<float>& dispR)
{
	// This function perform several times of weighted median filtering at each invalid position.
	// weights of the neighborhood are set by exp(-||cp-cq|| / gamma), except in the last iteration,
	// where the weights of invalid pixels are set to zero.

	VECBITMAP<bool> validL(nrows, ncols), validR(nrows, ncols);
	PlaneMapToDisparityMap(coeffsL, dispL);
	PlaneMapToDisparityMap(coeffsR, dispR);

#ifdef DO_POST_PROCESSING
	// Hole filling
	CrossCheck(dispL, dispR, validL, validR);
	#pragma omp parallel for
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			if (!validL[y][x]) {
				FillHole(y, x, validL, coeffsL);
			}
			if (!validR[y][x]) {
				FillHole(y, x, validR, coeffsR);
			}
		}
	}

	// Weighted median filtering 
	int maxround = 1;
	bool useInvalidPixels = true;
	for (int round = 0; round < maxround; round++) {

		PlaneMapToDisparityMap(coeffsL, dispL);
		PlaneMapToDisparityMap(coeffsR, dispR);
		CrossCheck(dispL, dispR, validL, validR);

		//if (round + 1 == maxround) {
		//	useInvalidPixels = false;
		//}
		#pragma omp parallel for
		for (int y = 0; y < nrows; y++) {
			//if (y % 10 == 0) { printf("median filtering row %d\n", y); }
			for (int x = 0; x < ncols; x++) {
				if (!validL[y][x]){
					VECBITMAP<float> wL(patch_w, patch_w, 1, weightsL.line_n1(y * ncols + x));
					WeightedMedianFilter(y, x, dispL, wL, validL, useInvalidPixels);
				}
				if (!validR[y][x]){
					VECBITMAP<float> wR(patch_w, patch_w, 1, weightsR.line_n1(y * ncols + x));
					WeightedMedianFilter(y, x, dispR, wR, validR, useInvalidPixels);
				}
			}
		}
	}
#endif
}

void RunPatchMatchStereo(cv::Mat& imL, cv::Mat& imR, int ndisps)
{
	VECBITMAP<float> dsiL = ComputeAdGradientCostVolume(imL, imR, ndisps, -1, granularity);
	VECBITMAP<float> dsiR = ComputeAdGradientCostVolume(imR, imL, ndisps, +1, granularity);
	VECBITMAP<float> weightsL = PrecomputeWeights(imL);
	VECBITMAP<float> weightsR = PrecomputeWeights(imR);

	VECBITMAP<float> dispL(nrows, ncols), dispR(nrows, ncols);
	VECBITMAP<Plane> coeffsL(nrows, ncols), coeffsR(nrows, ncols);
	VECBITMAP<float> bestcostsL(nrows, ncols), bestcostsR(nrows, ncols);

#ifndef LOAD_RESULT_FROM_LAST_RUN
	// Random initialization
	Timer::tic("Random Init");
	RandomInit(coeffsL, bestcostsL, dsiL, weightsL);
	RandomInit(coeffsR, bestcostsR, dsiR, weightsR);
	Timer::toc();

	// Iteration
	for (int iter = 0; iter < maxiters; iter++) {

		if (iter % 2 == 0) {
			Timer::tic("Left View");
			#pragma omp parallel for
			for (int y = 0; y < nrows; y++) {
				for (int x = 0; x < ncols; x++) {
					PropagateAndRandomSearch(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, dsiL, dsiR, weightsL, weightsR, iter, -1);
				}
			}
			Timer::toc();
			Timer::tic("Right View");
			#pragma omp parallel for
			for (int y = 0; y < nrows; y++) {
				for (int x = 0; x < ncols; x++) {
					PropagateAndRandomSearch(y, x, coeffsR, coeffsL, bestcostsR, bestcostsL, dsiR, dsiL, weightsR, weightsL, iter, +1);
				}
			}
			Timer::toc();
		}
		else {
			Timer::tic("Left View");
			#pragma omp parallel for
			for (int y = nrows - 1; y >= 0; y--) {
				for (int x = ncols - 1; x >= 0; x--) {
					PropagateAndRandomSearch(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, dsiL, dsiR, weightsL, weightsR, iter, -1);
				}
			}
			Timer::toc();
			Timer::tic("Right View");
			#pragma omp parallel for
			for (int y = nrows - 1; y >= 0; y--) {
				for (int x = ncols - 1; x >= 0; x--) {
					PropagateAndRandomSearch(y, x, coeffsR, coeffsL, bestcostsR, bestcostsL, dsiR, dsiL, weightsR, weightsL, iter, +1);
				}
			}
			Timer::toc();
		}
	}

	printf("g_improve_cnt: %d\n", g_improve_cnt);
	coeffsL.SaveToBinaryFile(folders[folder_id] + "coeffsL.bin");
	coeffsR.SaveToBinaryFile(folders[folder_id] + "coeffsR.bin");
#else
	coeffsL.LoadFromBinaryFile(folders[folder_id] + "coeffsL.bin");
	coeffsR.LoadFromBinaryFile(folders[folder_id] + "coeffsR.bin");
#endif

	// Post processing
	Timer::tic("POST PROCESSING");
	PostProcess(weightsL, weightsR, coeffsL, coeffsR, dispL, dispR);
	Timer::toc();

	EvaluateDisparity(dispL, 0.5f, coeffsL);
}




int main()
{
#ifndef USE_OPENMP
	omp_set_num_threads(1);
#endif

	cv::Mat imL = cv::imread(folders[folder_id] + "im2.png");
	cv::Mat imR = cv::imread(folders[folder_id] + "im6.png");
	nrows = imL.rows;
	ncols = imL.cols;


	cv::Mat dispL, dispR;
	//Timer::tic("LocalSearch");
	//LocalSearch(imL, imR, ndisps, dispL, dispR);
	Timer::tic("PatchMatchStereo");
	RunPatchMatchStereo(imL, imR, ndisps);
	//RunLaplacianStereo(imL, imR, ndisps);
	Timer::toc();

	
	//dispL.convertTo(dispL, CV_8UC3, scale);
	//dispR.convertTo(dispL, CV_8UC3, scale);
	//cv::imshow("disparity", dispL);
	//cv::waitKey(0);

	return 0;
}

