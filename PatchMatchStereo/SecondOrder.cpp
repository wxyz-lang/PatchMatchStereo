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

#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <omp.h>
#include "Utilities.h"


extern int nrows, ncols;


Eigen::SparseMatrix<double> PrecomputeSparseLTL(cv::Mat& cvImg)
{
	assert(cvImg.isContinuous());
	VECBITMAP<unsigned char> img(nrows, ncols, 3, cvImg.data);
	cv::Mat weightImgLR(nrows, ncols, CV_32FC1);
	cv::Mat weightImgUD(nrows, ncols, CV_32FC1);

	// Construct LTWT.
	const double sigma = 90.f;
	const int N = nrows * ncols;
	std::vector<Eigen::Triplet<double>> coefficientsLR;
	std::vector<Eigen::Triplet<double>> coefficientsUD;
	coefficientsLR.reserve(3 * N);
	coefficientsUD.reserve(3 * N);

	for (int y = 0, i = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++, i++) {
			if (x > 0 && x < ncols - 1) {
				Eigen::Triplet<double> left(i, i - 1, -1);
				Eigen::Triplet<double> right(i, i + 1, -1);
				Eigen::Triplet<double> self(i, i, 2);
				coefficientsLR.push_back(left);
				coefficientsLR.push_back(right);
				coefficientsLR.push_back(self);
			}
			if (y > 0 && y < nrows - 1) {
				Eigen::Triplet<double> up(i, i - ncols, -1);
				Eigen::Triplet<double> down(i, i + ncols, -1);
				Eigen::Triplet<double> self(i, i, 2);
				coefficientsUD.push_back(up);
				coefficientsUD.push_back(down);
				coefficientsUD.push_back(self);
			}
		}
	}

	Eigen::SparseMatrix<double> L1(N, N);
	Eigen::SparseMatrix<double> L2(N, N);
	L1.setFromTriplets(coefficientsLR.begin(), coefficientsLR.end());
	L2.setFromTriplets(coefficientsUD.begin(), coefficientsUD.end());

	coefficientsLR.clear();
	coefficientsUD.clear();
	coefficientsLR.reserve(1 * N);
	coefficientsUD.reserve(1 * N);
	for (int y = 0, i = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++, i++) {
			if (x > 0 && x < ncols - 1) {
				unsigned char *p = img.get(y, x);
				unsigned char *pL = img.get(y, x - 1);
				unsigned char *pR = img.get(y, x + 1);
				float diffL = fabs((float)p[0] - pL[0]) + fabs((float)p[1] - pL[1]) + fabs((float)p[2] - pL[2]);
				float diffR = fabs((float)p[0] - pR[0]) + fabs((float)p[1] - pR[1]) + fabs((float)p[2] - pR[2]);
				float wL = exp(-diffL / sigma);
				float wR = exp(-diffR / sigma);
				float weight = std::min(wL, wR);

				Eigen::Triplet<double> triplet(i, i, weight);
				coefficientsLR.push_back(triplet);

				weightImgLR.at<float>(y, x) = weight;
			}
			if (y > 0 && y < nrows - 1) {
				unsigned char *p = img.get(y, x);
				unsigned char *pU = img.get(y - 1, x);
				unsigned char *pD = img.get(y + 1, x);
				float diffU = fabs((float)p[0] - pU[0]) + fabs((float)p[1] - pU[1]) + fabs((float)p[2] - pU[2]);
				float diffD = fabs((float)p[0] - pD[0]) + fabs((float)p[1] - pD[1]) + fabs((float)p[2] - pD[2]);
				float wU = exp(-diffU / sigma);
				float wD = exp(-diffD / sigma);
				float weight = std::min(wU, wD);

				Eigen::Triplet<double> triplet(i, i, weight);
				coefficientsUD.push_back(triplet);

				weightImgUD.at<float>(y, x) = weight;
			}
		}
	}

	//weightImgLR.convertTo(weightImgLR, CV_8UC3, 255);
	//weightImgUD.convertTo(weightImgUD, CV_8UC3, 255);
	//cv::imshow("weightLR", weightImgLR);
	//cv::imshow("weightUD", weightImgUD);
	//cv::waitKey(0);

	Eigen::SparseMatrix<double> W1(N, N);
	Eigen::SparseMatrix<double> W2(N, N);
	W1.setFromTriplets(coefficientsLR.begin(), coefficientsLR.end());
	W2.setFromTriplets(coefficientsUD.begin(), coefficientsUD.end());

	Eigen::SparseMatrix<double> LTL = L1.transpose() * W1 * L1 + L2.transpose() * W2 * L2;
	printf("L.nonZeros = %d, LTL.nonZeros = %d\n", L1.nonZeros(), LTL.nonZeros());

	return LTL;
}

VECBITMAP<float> SolveSecondOrderSmootheness(VECBITMAP<float>& dv, float theta, Eigen::SparseMatrix<double>& LTL)
{
	const int N = nrows * ncols;
	Eigen::SparseMatrix<double> G(N, N);
	G.setIdentity();

	Eigen::VectorXd v(N);
	for (int i = 0; i < N; i++) {
		v[i] = dv.data[i];
	}

	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(LTL + theta * G);
	Eigen::VectorXd u = chol.solve(theta * G * v);
	
	VECBITMAP<float> du(nrows, ncols);
	for (int i = 0; i < N; i++) {
		du.data[i] = u[i];
	}

	return du;
}

VECBITMAP<float> ConstrainedLocalSearch(VECBITMAP<float>&u, VECBITMAP<float>& dsi, float theta, float lambda)
{
	VECBITMAP<float> dsi_constrained(nrows, ncols, ndisps);
	memcpy(dsi_constrained.data, dsi.data, nrows * ncols * ndisps * sizeof(float));

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			for (int d = 0; d < ndisps; d++) {
				dsi_constrained.get(y, x)[d] *= lambda;
				dsi_constrained.get(y, x)[d] += theta * (d - u[y][x]) * (d - u[y][x]);
			}
		}
	}

	return WinnerTakesAll(dsi_constrained);
}

void RunLaplacianStereo(cv::Mat& imL, cv::Mat& imR, int ndisps)
{
	// Initialize u and v from GT.
	cv::Mat gt = cv::imread(folders[folder_id] + "disp2.png", CV_LOAD_IMAGE_GRAYSCALE);
	gt.convertTo(gt, CV_32FC1);
	gt /= 4;
	assert(gt.isContinuous());

	//VECBITMAP<float> u(nrows, ncols), v(nrows, ncols);
	//memcpy(u.data, gt.data, nrows * ncols * sizeof(float));
	//memcpy(v.data, gt.data, nrows * ncols * sizeof(float));

	//VECBITMAP<float> dsiL = ComputeAdGradientCostVolume(imL, imR, ndisps, -1, 1.f);
	VECBITMAP<float> dsiL = ComputeAdCensusCostVolume(imL, imR, ndisps, -1);
	VECBITMAP<float> u = WinnerTakesAll(dsiL);
	VECBITMAP<float> v = WinnerTakesAll(dsiL);

	Timer::tic("Prepare LTL matrix");
	Eigen::SparseMatrix<double> LTL = PrecomputeSparseLTL(imL);
	Timer::toc();

	
	const float lambda = 20.f;

	EvaluateDisparity(u, 0.5f);

	// Iterate between smoothness and matching cost
	for (float theta = 0.001; theta < 20; theta *= 1.5f) {

		printf("theta = %f\n", theta);

		Timer::tic("SolveSecondOrderSmoothness");
		u = SolveSecondOrderSmootheness(v, theta, LTL);
		Timer::toc();
		//EvaluateDisparity(u, 0.5f);

		Timer::tic("ConstrainedLocalSearch");
		v = ConstrainedLocalSearch(u, dsiL, theta, lambda);
		Timer::toc();
		//EvaluateDisparity(v, 0.5f);
	}

	//for (float theta = 0.001; theta < 20; theta *= 1.5f) {

	//	printf("theta = %f\n", theta);

	//	Timer::tic("SolveSecondOrderSmoothness");
	//	uL = SolveSecondOrderSmootheness(vL, theta, LTL);
	//	uR = SolveSecondOrderSmootheness(vR, theta, LTL);
	//	Timer::toc();
	//	//EvaluateDisparity(u, 0.5f);

	//	Timer::tic("ConstrainedLocalSearch");
	//	RunPatchMatchStereo(imL, imR, ndisps, uL, uR);
	//	Timer::toc();
	//	//EvaluateDisparity(v, 0.5f);
	//}

	EvaluateDisparity(u, 0.5f);
	EvaluateDisparity(v, 0.5f);
}
