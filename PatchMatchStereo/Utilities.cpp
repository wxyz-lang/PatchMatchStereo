
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





extern int nrows, ncols;

cv::Mat g_L, g_R, g_OC, g_ALL, g_DISC, g_GT, g_segments;
cv::Mat g_stereo, g_mydisp;
VECBITMAP<Plane> g_coeffsL;
VECBITMAP<float> g_colgradL, g_colgradR;

extern std::vector<std::vector<cv::Point2d>> g_regionList;
extern VECBITMAP<int> g_labelmap;
extern VECBITMAP<Plane> g_coeffsL_ransac, g_coeffsL_neldermead;
extern VECBITMAP<float> g_dsiL;


VECBITMAP<float> ComputeColGradFeature(cv::Mat& img);



inline bool InBound(float y, float x) { return 0 <= y && y < nrows && 0 <= x && x < ncols; }

void prepareRect(int yc, int xc, cv::Rect& patchRect, cv::Rect& canvasRect)
{
	int yU = std::max(0, yc - patch_r);
	int xL = std::max(0, xc - patch_r);
	int yD = std::min(nrows - 1, yc + patch_r);
	int xR = std::min(ncols - 1, xc + patch_r);
	canvasRect = cv::Rect(xL, yU, xR - xL + 1, yD - yU + 1);
	// From canvas coordinate to patch coordinate
	patchRect = cv::Rect(xL - xc + patch_r, yU - yc + patch_r, xR - xL + 1, yD - yU + 1);
}

VECBITMAP<float> ComputeLocalPatchWeights(int yc, int xc, cv::Mat& cvImg)
{
	assert(cvImg.isContinuous());
	VECBITMAP<unsigned char> img(nrows, ncols, 3, cvImg.data);
	VECBITMAP<float> ret(patch_w, patch_w);
	float *w = ret.data;
	memset(w, 0, patch_w * patch_w * sizeof(float));
	unsigned char *rgb1 = img.get(yc, xc);
	unsigned char *rgb2 = NULL;

	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

	#pragma omp parallel for
	for (int y = yb; y <= ye; y++) {
		for (int x = xb; x <= xe; x++) {
			rgb2 = img.get(y, x);
			float cost = std::abs((float)rgb1[0] - (float)rgb2[0])
				+ std::abs((float)rgb1[1] - (float)rgb2[1])
				+ std::abs((float)rgb1[2] - (float)rgb2[2]);
			w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] = cost;
		}
	}
	for (int i = 0; i < patch_w * patch_w; i++) {
		w[i] = exp(-w[i] / gamma);
	}
	return ret;
}

double ComputePlaneCost(int yc, int xc, Plane& coeff_try, VECBITMAP<float>& colgradL, VECBITMAP<float> &colgradR, VECBITMAP<float>& w, int sign)
{
	const float BAD_PLANE_PENALTY = 3 * ((1 - alpha) * tau_col + alpha * tau_grad);		// set to 3 times the max possible matchin cost.
	double cost = 0;
	for (int y = yc - patch_r; y <= yc + patch_r; y++) {
		for (int x = xc - patch_r; x <= xc + patch_r; x++) {
			//int d = 0.5 + (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
			float d = (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
			if (InBound(y, x)) {


				if (d < 0 || d > dmax) {	// must be a bad plane.
					cost += BAD_PLANE_PENALTY;
				}
				else {

					float xm, xmL, xmR, wL, wR;
					float cost_col, cost_grad;
					xm = std::max(0.f, std::min((float)ncols - 1, x + sign * d));
					xmL = (int)(xm);
					xmR = (int)(xm + 0.99);
					wL = xmR - xm;
					wR = 1 - wL;

					float *pL = colgradL.get(y, x);
					float *pRmL = colgradR.get(y, xmL);
					float *pRmR = colgradR.get(y, xmR);
					float pR[5];
					for (int i = 0; i < 5; i++) {
						pR[i] = wL * pRmL[i] + wR * pRmR[i];
					}
					//float *pR = colgradR.get(y, (int)(xm + 0.5));
					cost_col = fabs(pL[0] - pR[0])
						+ fabs(pL[1] - pR[1])
						+ fabs(pL[2] - pR[2]);
					cost_col = std::min(tau_col, cost_col);
					cost_grad = fabs(pL[3] - pR[3])
						+ fabs(pL[4] - pR[4]);
					cost_grad = std::min(tau_grad, cost_grad);
					cost += w[y - yc + patch_r][x - xc + patch_r] * ((1 - alpha) * cost_col + alpha * cost_grad);
				}
			}
		}
	}
	return cost;
}

void on_mouse(int event, int x, int y, int flags, void *param)
{

	CvPoint left, right;
	cv::Mat tmp;

	if (event == CV_EVENT_MOUSEMOVE)
	{
		if (x >= ncols && x < 2 * ncols && y >= 0 && y < 2 * nrows)
		{
			tmp = g_stereo.clone();
			left.x = x; left.y = y;
			right.x = x + ncols; 	right.y = y;
			cv::line(tmp, left, right, cv::Scalar(255, 0, 0, 1));

			left.x = x; left.y = y;
			right.x = x - ncols; 	right.y = y;
			cv::line(tmp, left, right, cv::Scalar(255, 0, 0, 1));

			cv::imshow("disparity", tmp);
		}
	}
	if (event == CV_EVENT_RBUTTONDOWN)
	{
		if (x >= ncols && x < 3 * ncols)
		{
			tmp = g_stereo.clone();
			x %= ncols; y %= nrows;
			float gtDisp = g_GT.at<cv::Vec3b>(y, x)[0];
			float myDisp = g_mydisp.at<cv::Vec3b>(y, x)[0];
			printf("(%d, %d)->GT: %.1f, MINE: %.1f    ", y, x, gtDisp / scale, myDisp / scale);
			if (g_coeffsL.data) {
				printf("a = %.5f  b = %.5f  c = %.5f\n", g_coeffsL[y][x].a, g_coeffsL[y][x].b, g_coeffsL[y][x].c);
			}

			VECBITMAP<float> weights = ComputeLocalPatchWeights(y, x, g_L);
			float wsum = 0;  for (int i = 0; i < patch_w * patch_w; i++) wsum += weights.data[i];
			if (g_coeffsL.data) {
				float cost = ComputePlaneCost(y, x, g_coeffsL[y][x], g_colgradL, g_colgradR, weights, -1) / wsum;
				printf("plane cost: %f\n", cost);
			}
			if (g_regionList.size() > 0) {

				int id = g_labelmap[y][x];
				Plane coeff_ransac = g_coeffsL_ransac[y][x];
				Plane coeff_neldermead = g_coeffsL_neldermead[y][x];

				double ComputePlaneCost(Plane& coeff, VECBITMAP<float>& dsi, std::vector<cv::Point2d>& pointList);
				float cost_ransac = ComputePlaneCost(coeff_ransac, g_dsiL, g_regionList[id]);
				float cost_neldermead = ComputePlaneCost(coeff_neldermead, g_dsiL, g_regionList[id]);

				printf("     RANSAC plane cost: %f\n", cost_ransac);
				printf("NELDER-MEAD plane cost: %f\n", cost_neldermead);
			}
			

			cv::Mat outImg1 = tmp(cv::Rect(0, 0, ncols, nrows)),
				outImg2 = tmp(cv::Rect(ncols, 0, ncols, nrows)),
				outImg3 = tmp(cv::Rect(0, nrows, ncols, nrows)),
				outImg4 = tmp(cv::Rect(ncols, nrows, ncols, nrows)),
				outImg5 = tmp(cv::Rect(ncols * 2, 0, ncols, nrows)),
				outImg6 = tmp(cv::Rect(ncols * 2, nrows, ncols, nrows));
			cv::Rect patchRect, canvasRect;
			prepareRect(y, x, patchRect, canvasRect);

			cv::Mat weightImg(patch_w, patch_w, CV_8UC3);
			for (int y = 0; y < patch_w; y++) {
				for (int x = 0; x < patch_w; x++) {
					unsigned char w = weights[y][x] * 255.f;
					weightImg.at<cv::Vec3b>(y, x) = cv::Vec3b(w, w, w);
				}
			}
			(weightImg(patchRect)).copyTo(outImg5(canvasRect));

			Plane coeff = g_coeffsL[y][x];
			int yc = y, xc = x;
			cv::Mat planeImg(patch_w, patch_w, CV_8UC3);
			for (int y = 0; y < patch_w; y++) {
				for (int x = 0; x < patch_w; x++) {
					unsigned char d = scale * coeff.ToDisparity(y - patch_r + yc, x - patch_r + xc);
					planeImg.at<cv::Vec3b>(y, x) = cv::Vec3b(d, d, d);
				}
			}
			(planeImg(patchRect)).copyTo(outImg2(canvasRect));
			cv::imshow("disparity", tmp);
		}
	}
}

void EvaluateDisparity(VECBITMAP<float>& h_disp, float thresh, VECBITMAP<Plane>& coeffsL)
{
	g_L = cv::imread(folders[folder_id] + "im2.png");
	g_R = cv::imread(folders[folder_id] + "im6.png");
	g_GT = cv::imread(folders[folder_id] + "disp2.png");
	g_OC = cv::imread(folders[folder_id] + "nonocc.png");
	g_ALL = cv::imread(folders[folder_id] + "all.png");
	g_DISC = cv::imread(folders[folder_id] + "disc.png");

	float count[3] = { 0, 0, 0 };
	float badPixelRate[3] = { 0, 0, 0 };
	cv::Mat disp(nrows, ncols, CV_8UC3), badOnALL(nrows, ncols, CV_8UC3), badOnOC(nrows, ncols, CV_8UC3), gray;
	cv::cvtColor(g_L, gray, CV_BGR2GRAY);

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {

			count[0] += (g_OC.at<cv::Vec3b>(y, x)[0] == 255);
			count[1] += (g_ALL.at<cv::Vec3b>(y, x)[0] == 255);
			count[2] += (g_DISC.at<cv::Vec3b>(y, x)[0] == 255);

			unsigned char g = gray.at<unsigned char>(y, x);
			badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(g, g, g);
			badOnOC.at<cv::Vec3b>(y, x) = cv::Vec3b(g, g, g);

			unsigned char d = (unsigned int)(scale * h_disp[y][x] + 0.5);
			disp.at<cv::Vec3b>(y, x) = cv::Vec3b(d, d, d);

			float diff = abs((float)disp.at<cv::Vec3b>(y, x)[0] - (float)g_GT.at<cv::Vec3b>(y, x)[0]);
			if (g_OC.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh) {
				badPixelRate[0]++;
				badOnOC.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 255);
			}
			if (g_ALL.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh) {
				badPixelRate[1]++;
				badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 255);
			}
			if (g_DISC.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh)	{
				badPixelRate[2]++;
			}



			if (g_OC.at<cv::Vec3b>(y, x)[0] == 0) {		//draw occlusion region
				badOnOC.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 0, 0);
			}
			if (g_ALL.at<cv::Vec3b>(y, x)[0] == 0) {	//draw occlusion region
				badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 0, 0);
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		badPixelRate[i] /= count[i];
	}

	printf("badPixelRate: %f%%, %f%%, %f%%\n",
		badPixelRate[0] * 100.0f, badPixelRate[1] * 100.0f, badPixelRate[2] * 100.0f);

	if (1)
	{

		cv::Mat compareImg(nrows * 2, ncols * 3, CV_8UC3);
		cv::Mat outImg1 = compareImg(cv::Rect(0, 0, ncols, nrows)),
			outImg2 = compareImg(cv::Rect(ncols, 0, ncols, nrows)),
			outImg3 = compareImg(cv::Rect(0, nrows, ncols, nrows)),
			outImg4 = compareImg(cv::Rect(ncols, nrows, ncols, nrows)),
			outImg5 = compareImg(cv::Rect(ncols * 2, 0, ncols, nrows)),
			outImg6 = compareImg(cv::Rect(ncols * 2, nrows, ncols, nrows));
		g_GT.copyTo(outImg1);	disp.copyTo(outImg2);
		g_segments.copyTo(outImg3);	badOnOC.copyTo(outImg4);
		g_L.copyTo(outImg5);	badOnALL.copyTo(outImg6);


		std::string folderpath = folders[folder_id];
		char buf[256];
		//sprintf(buf, "%s=%.2f.png", folderpath.substr(0, folderpath.length() - 1).c_str(), badPixelRate[0] * 100.0f);
		//std::string filename(buf);
		std::string filename = folderpath.substr(0, folderpath.length() - 1) + ".png";
		std::string filename2 = folderpath.substr(0, folderpath.length() - 1) + "_err.png";
		cv::imwrite("d:\\" + filename, disp);
		cv::imwrite("d:\\" + filename2, compareImg);

		cv::imwrite(folders[folder_id] + "disparity.png", disp);
		cv::imwrite("d:\\disparity.png", disp);
		cv::imwrite("d:\\badOnALL.png", badOnALL);
		cv::imwrite("d:\\badOnOC.png", badOnOC);

		g_colgradL = ComputeColGradFeature(g_L);
		g_colgradR = ComputeColGradFeature(g_R);
		g_coeffsL = coeffsL;
		g_stereo = compareImg;
		g_mydisp = disp;

		cv::imshow("disparity", compareImg);
		cv::setMouseCallback("disparity", on_mouse);
		cv::waitKey(0);
	}
}