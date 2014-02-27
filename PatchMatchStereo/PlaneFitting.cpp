
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
#include "msImageProcessor.h"
#include "SLIC.h"


extern int nrows, ncols;
std::vector<std::vector<cv::Point2d>> g_regionList;
VECBITMAP<int> g_labelmap;
VECBITMAP<Plane> g_coeffsL_ransac, g_coeffsL_neldermead;
VECBITMAP<float> g_dsiL;


int meanShiftSegmentation(const cv::Mat &img, const float colorRadius, const int spatialRadius, const int minRegion, cv::Mat& result)
{
	const int width = img.cols, height = img.rows;
	cv::Mat rgb;
	cv::cvtColor(img, rgb, cv::COLOR_BGR2RGB);

	msImageProcessor ms;
	ms.DefineImage((unsigned char*)rgb.data, COLOR, height, width);
	if (ms.ErrorStatus)
	{
		return -1;
	}

	ms.Filter(spatialRadius, colorRadius, NO_SPEEDUP);
	ms.FuseRegions(colorRadius, minRegion);
	if (ms.ErrorStatus)
	{
		return -1;
	}

	cv::Mat segmImage(height, width, CV_8UC3);
	ms.GetResults(segmImage.data);
	if (ms.ErrorStatus)
	{
		return -1;
	}
	cv::cvtColor(segmImage, result, CV_RGB2BGR);
	//cv::imshow("segments", result);
	//cv::waitKey(0);

	const cv::Scalar colorDiff = cv::Scalar::all(0);
	cv::Mat mask(result.rows + 2, result.cols + 2, CV_8UC1, cv::Scalar::all(0));
	for (int y = 0; y < result.rows; y++)
	{
		for (int x = 0; x < result.cols; x++)
		{
			if (mask.at<uchar>(y, x) == 0)
			{
				cv::Scalar newVal(rand() % 255, rand() % 255, rand() % 255);
				cv::floodFill(result, mask, cv::Point(x, y), newVal, 0, colorDiff, colorDiff);
			}
		}
	}
	//cv::imwrite("d:\\segments.png", result);
	//cv::imshow("meanshift", result);
	//cv::waitKey(0);

	return 1;
}

int slicSegmentation(const cv::Mat &img, const int numPreferedRegions, const int compactness, cv::Mat& result)
{
	cv::Mat argb(nrows, ncols, CV_8UC4);
	assert(argb.isContinuous());

	int from_to[] = { -1, 0, 0, 3, 1, 2, 2, 1 };
	cv::mixChannels(&img, 1, &argb, 1, from_to, 4);

	int width(ncols), height(nrows), numlabels(0);;
	unsigned int* pbuff = (unsigned int*)argb.data;
	int* klabels = NULL;

	int		k = numPreferedRegions;	// Desired number of superpixels.
	double	m = compactness;		// Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
	SLIC segment;
	segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, width, height, klabels, numlabels, k, m);
	//segment.DrawContoursAroundSegments(pbuff, klabels, width, height, 0xff0000);
	VECBITMAP<int> labelmap(nrows, ncols, 1, klabels);

	//printf("numlabels = %d\n", numlabels);
	
	std::vector<cv::Vec3b> colors(numlabels);
	for (int i = 0; i < numlabels; i++) {
		colors[i] = cv::Vec3b(rand() % 255, rand() % 255, rand() % 255);
	}
	result.create(nrows, ncols, CV_8UC3);
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			int id = labelmap[y][x];
			result.at<cv::Vec3b>(y, x) = colors[id];
		}
	}
	
	delete[] klabels;

	//cv::imwrite("d:\\segments.png", result);
	//cv::imshow("meanshift", result);
	//cv::waitKey(0);

	return 1;
}

inline bool InBound(float y, float x) { return 0 <= y && y < nrows && 0 <= x && x < ncols; }

int FromSegmentMapToLabelMap(cv::Mat& segmap, VECBITMAP<int>& labelmap)
{
	const int dx[] = { -1, 0, +1, 0 };
	const int dy[] = { 0, -1, 0, +1 };

	VECBITMAP<bool> visited(nrows, ncols);
	memset(visited.data, 0, nrows * ncols * sizeof(bool));
	int label = -1;
	std::stack<cv::Point2d> stack;

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			if (!visited[y][x]) {
				label++;
				stack.push(cv::Point2d(x, y));
				
				while (!stack.empty()) {
					cv::Point2d c = stack.top();
					stack.pop();
					int yc = c.y;
					int xc = c.x;
					labelmap[yc][xc] = label;
					visited[yc][xc] = true;
					for (int dir = 0; dir < 4; dir++) {
						int new_x = xc + dx[dir];
						int new_y = yc + dy[dir];
						if (InBound(new_y, new_x) && !visited[new_y][new_x] && segmap.at<cv::Vec3b>(new_y, new_x) == segmap.at<cv::Vec3b>(yc, xc)) {
							stack.push(cv::Point2d(new_x, new_y));
						}
					}
				}
			}
		}
	}

	return label + 1;
}

void RandomPermute(std::vector<cv::Point2d>& pointList, int k)
{
	int n = pointList.size();
	for (int i = 0; i < k; i++) {
		int j = rand() % n;
		std::swap(pointList[i], pointList[j]);
	}
}

double ComputePlaneCost(Plane& coeff, VECBITMAP<float>& dsi, std::vector<cv::Point2d>& pointList)
{
	double cost = 0;
	int regionSize = pointList.size();
	for (int i = 0; i < regionSize; i++) {
		int y = pointList[i].y;
		int x = pointList[i].x;
		float d = coeff.a * x + coeff.b * y + coeff.c;
		int level = d / granularity + 0.5;
		if (d < 0 || d > dmax) {
			cost += BAD_PLANE_PENALTY;
		}
		else {
			cost += dsi.get(y, x)[level];
		}
	}
	return cost;
}

struct NM_OPT_PARAM {
	VECBITMAP<float> *dsi;
	std::vector<cv::Point2d> *pointList;
};
NM_OPT_PARAM nm_opt_struct;

inline float nm_compute_plane_cost(float *abc, int n = 3)
{
	Plane coeff;
	coeff.SetAbc(abc);
	return ComputePlaneCost(coeff, *nm_opt_struct.dsi, *nm_opt_struct.pointList);
}

void NelderMeadEstimate(std::vector<cv::Point2d>& pointList, VECBITMAP<float>& dsi, VECBITMAP<float>& disp, VECBITMAP<Plane>& coeffs)
{

	int NelderMeadOptimize(float *x, float(*feval)(float*, int), int maxiters = 0);

	const int MIN_SAMPLE_SIZE = 5;
	const int regionSize = pointList.size();

	if (regionSize < MIN_SAMPLE_SIZE) {
		for (int i = 0; i < regionSize; i++) {
			int y = pointList[i].y;
			int x = pointList[i].x;
			coeffs[y][x].RandomAssign(y, x, dmax);
		}
		return;
	}

	for (int i = 0; i < pointList.size(); i++) {
		int y = pointList[i].y;
		int x = pointList[i].x;
		float wta_d = disp[y][x];
		coeffs[y][x].RandomAssignNormal(y, x, dmax, wta_d);
	}


	float vertices[3 * 4];
	int max_iters = std::min(regionSize, 20);
	for (int iter = 0; iter < max_iters; iter++) {
		// initialize simplex vertices
		RandomPermute(pointList, 4);
		int i;
		if (iter == 0) { i = 0; }
		else { i = 1; }
		for (; i < 4; i++) {
			int y = pointList[i].y;
			int x = pointList[i].x;
			vertices[3 * i + 0] = coeffs[y][x].a;
			vertices[3 * i + 1] = coeffs[y][x].b;
			vertices[3 * i + 2] = coeffs[y][x].c;
		}

		//int y = pointList[0].y;
		//int x = pointList[0].x;
		//Plane coeff = coeffs[y][x];
		//vertices[0] = coeff.a;
		//vertices[1] = coeff.b;
		//vertices[2] = coeff.c;
		//for (int i = 1; i < 4; i++) {
		//	int y = pointList[i].y;
		//	int x = pointList[i].x;
		//	Plane coeff;
		//	coeff.RandomAssign(y, x, dmax);
		//	vertices[3 * i + 0] = coeff.a;
		//	vertices[3 * i + 1] = coeff.b;
		//	vertices[3 * i + 2] = coeff.c;
		//}

		// invoke nelder-mead
		nm_opt_struct.dsi = &dsi;
		nm_opt_struct.pointList = &pointList;
		float cost_before = nm_compute_plane_cost(vertices);
		NelderMeadOptimize(vertices, nm_compute_plane_cost, 30);
		float cost_after = nm_compute_plane_cost(vertices);

		if (cost_after - cost_before > 0) {
			printf("BUG: energy increased!\n");
		}
	}
	
	

	for (int i = 0; i < regionSize; i++) {
		int y = pointList[i].y;
		int x = pointList[i].x;
		coeffs[y][x].a = vertices[0];
		coeffs[y][x].b = vertices[1];
		coeffs[y][x].c = vertices[2];
	}
}

void RansacEstimate(std::vector<cv::Point2d>& pointList, VECBITMAP<float>& dsi, VECBITMAP<float>& disp, VECBITMAP<Plane>& coeffs)
{
	const int MIN_SAMPLE_SIZE = 5;
	const int regionSize = pointList.size();

	if (regionSize < MIN_SAMPLE_SIZE) {
		for (int i = 0; i < regionSize; i++) {
			int y = pointList[i].y;
			int x = pointList[i].x;
			coeffs[y][x].RandomAssign(y, x, dmax);
		}
		return;
	}

	cv::Mat A(MIN_SAMPLE_SIZE, 3, CV_64F);
	cv::Mat b(MIN_SAMPLE_SIZE, 1, CV_64F);
	cv::Mat c(3, 1, CV_64F);

	Plane bestcoeff;
	double bestcost = DBL_MAX;

	for (int retry = 0; retry < 200; retry++) {

		RandomPermute(pointList, MIN_SAMPLE_SIZE);
		for (int i = 0; i < MIN_SAMPLE_SIZE; i++) {
			int y = pointList[i].y;
			int x = pointList[i].x;
			A.at<double>(i, 0) = x;
			A.at<double>(i, 1) = y;
			A.at<double>(i, 2) = 1;
			b.at<double>(i, 0) = disp[y][x];
		}

		//printf("solving\n");
		cv::solve(A, b, c, cv::DECOMP_SVD);
		//printf("finished.\n");
		Plane coeff;
		coeff.a = c.at<double>(0, 0);
		coeff.b = c.at<double>(1, 0);
		coeff.c = c.at<double>(2, 0);

		double cost = ComputePlaneCost(coeff, dsi, pointList);
		if (cost < bestcost) {
			bestcost = cost;
			bestcoeff = coeff;
		}
	}

	for (int i = 0; i < regionSize; i++) {
		int y = pointList[i].y;
		int x = pointList[i].x;
		coeffs[y][x] = bestcoeff;
	}
}

void PlanefitView(cv::Mat& imL, VECBITMAP<float>& dsiL, VECBITMAP<Plane>& coeffsL, VECBITMAP<float>& dispL)
{
	extern cv::Mat g_segments;
	Timer::tic("segmentation");
	//meanShiftSegmentation(imL, 4, 5, 20, g_segments);
	slicSegmentation(imL, 300, 20, g_segments);
	Timer::toc();

	VECBITMAP<int> labelmap(nrows, ncols);
	int nlables = FromSegmentMapToLabelMap(g_segments, labelmap);


	std::vector<int> regionSize(nlables, 0);
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			regionSize[labelmap[y][x]]++;
		}
	}
	std::vector<std::vector<cv::Point2d>> regionList(nlables);
	for (int i = 0; i < nlables; i++) {
		regionList[i].reserve(regionSize[i]);
	}
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			int label = labelmap[y][x];
			regionList[label].push_back(cv::Point2d(x, y));
		}
	}

	g_regionList = regionList;
	g_labelmap = labelmap;
	g_dsiL = dsiL;

	Timer::tic("Fitting region");
	//#pragma omp parallel for
	for (int id = 0; id < nlables; id++) {
		//RansacEstimate(regionList[id], dsiL, dispL, coeffsL);
	}
	g_coeffsL_ransac = coeffsL;

	//#pragma omp parallel for
	for (int id = 0; id < nlables; id++) {	
		// This function is not thread-safe!
		NelderMeadEstimate(regionList[id], dsiL, dispL, coeffsL);
	}
	g_coeffsL_neldermead = coeffsL;
	Timer::toc();

	PlaneMapToDisparityMap(coeffsL, dispL);
}

void RunRansacPlaneFitting(cv::Mat& imL, cv::Mat& imR, int ndisps)
{
	VECBITMAP<float> dsiL = ComputeAdGradientCostVolume(imL, imR, ndisps, -1, granularity);
	VECBITMAP<float> dsiR = ComputeAdGradientCostVolume(imR, imL, ndisps, +1, granularity);

	VECBITMAP<float> adcensus_dsiL = ComputeAdCensusCostVolume(imL, imR, ndisps, -1);
	VECBITMAP<float> adcensus_dsiR = ComputeAdCensusCostVolume(imR, imL, ndisps, +1);

	//VECBITMAP<float> dispL = WinnerTakesAll(dsiL, granularity);
	//VECBITMAP<float> dispR = WinnerTakesAll(dsiR, granularity);
	VECBITMAP<float> dispL = WinnerTakesAll(adcensus_dsiL, 1.0);
	VECBITMAP<float> dispR = WinnerTakesAll(adcensus_dsiR, 1.0);
	//VECBITMAP<float> dispL(nrows, ncols), dispR(nrows, ncols);
	//cv::Mat gtL = cv::imread(folders[folder_id] + "disp2.png", CV_LOAD_IMAGE_GRAYSCALE);
	//cv::Mat gtR = cv::imread(folders[folder_id] + "disp6.png", CV_LOAD_IMAGE_GRAYSCALE);
	//gtL.convertTo(gtL, CV_32FC1, 0.25);
	//gtR.convertTo(gtR, CV_32FC1, 0.25);
	//assert(gtL.isContinuous());
	//assert(gtR.isContinuous());
	//memcpy(dispL.data, gtL.data, nrows*ncols*sizeof(float));
	//memcpy(dispR.data, gtR.data, nrows*ncols*sizeof(float));


	//cv::Mat dL(nrows, ncols, CV_32FC1);
	//memcpy(dL.data, dispL.data, nrows * ncols * sizeof(float));
	//dL.convertTo(dL, CV_8UC3, scale);
	//imshow("WTA", dL);
	//cv::waitKey(0);

	VECBITMAP<Plane> coeffsL(nrows, ncols);
	VECBITMAP<Plane> coeffsR(nrows, ncols);
	PlanefitView(imL, dsiL, coeffsL, dispL);
	PlanefitView(imR, dsiR, coeffsR, dispR);

	VECBITMAP<float> weightsL = PrecomputeWeights(imL);
	VECBITMAP<float> weightsR = PrecomputeWeights(imR);

	Timer::tic("Planefit Postprocess");
	PostProcess(weightsL, weightsR, coeffsL, coeffsR, dispL, dispR);
	Timer::toc();

	EvaluateDisparity(dispL, 0.5f, coeffsL);
}
