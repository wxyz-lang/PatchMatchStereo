
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

inline double ranf()	// ranf() is uniform in 0..1
{
	double a = rand();
	return a / RAND_MAX;
}

double BoxMuller(double m, double s)	// normal random variate generator 
{										// mean m, standard deviation s 
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)										// use value from previous call 
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	return(m + y1 * s);
}

struct QuadraticSurface {
	float A, B, C, D, E, F;
	QuadraticSurface() {}
	QuadraticSurface(float A_, float B_, float C_, float D_, float E_, float F_)
		: A(A_), B(B_), C(C_), D(D_), E(E_), F(F_) {}
	void RandomPerturbTo(float *buf)
	{
		buf[0] = BoxMuller(A, 0.001);
		buf[1] = BoxMuller(B, 0.001);
		buf[2] = BoxMuller(C, 0);
		buf[3] = BoxMuller(D, 0.01);
		buf[4] = BoxMuller(E, 0.01);
		buf[5] = BoxMuller(F, 3);
		//buf[0] = BoxMuller(A, 0.001);
		//buf[1] = BoxMuller(B, 0.001);
		//buf[2] = BoxMuller(C, 0.000);
		//buf[3] = BoxMuller(D, 0.00);
		//buf[4] = BoxMuller(E, 0.00);
		//buf[5] = BoxMuller(F, 0.00);
	}
	void SetFromArray(float *arr) {
		A = arr[0];
		B = arr[1];
		C = arr[2];
		D = arr[3];
		E = arr[4];
		F = arr[5];
	}
	void SetToArray(float *arr)
	{
		arr[0] = A;
		arr[1] = B;
		arr[2] = C;
		arr[3] = D;
		arr[4] = E;
		arr[5] = F;
	}
	float ToDisparity(float y, float x)
	{
		return A*x*x + B*y*y + C*x*y + D*x + E*y + F;
	}
};

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
	float x0, y0;
	VECBITMAP<float> *dsi;
	std::vector<cv::Point2d> *pointList;
};
NM_OPT_PARAM nm_opt_struct;


double ComputeQuadraticSurfaceCost(QuadraticSurface& coeff, VECBITMAP<float>& dsi, std::vector<cv::Point2d>& pointList)
{
	double cost = 0;
	int regionSize = pointList.size();
	int x0 = nm_opt_struct.x0;
	int y0 = nm_opt_struct.y0;
	//printf("regionSize = %d\n", regionSize);
	for (int i = 0; i < regionSize; i++) {
		int y = pointList[i].y;
		int x = pointList[i].x;
		float d = coeff.A * x * x + coeff.B * y * y + coeff.C * x * y
		        + coeff.D * x     + coeff.E * y     + coeff.F;
		int level = d / granularity + 0.5;
		//printf("(%d, %d)   i = %d   d = %f   level = %d\n", y, x, i, d, level);
		if (d < 0 || d > dmax) {
			cost += BAD_PLANE_PENALTY;
		}
		else {
			cost += dsi.get(y + y0, x + x0)[level];
		}
	}
	return cost;
}


inline float nm_compute_plane_cost(float *abc, int n = 3)
{
	Plane coeff;
	coeff.SetAbc(abc);
	return ComputePlaneCost(coeff, *nm_opt_struct.dsi, *nm_opt_struct.pointList);
}

inline float nm_compute_quadratic_surface_cost(float *abcdef, int n = 6)
{
	QuadraticSurface coeff;
	coeff.SetFromArray(abcdef);
	//printf("ready\n");
	return ComputeQuadraticSurfaceCost(coeff, *nm_opt_struct.dsi, *nm_opt_struct.pointList);
}


void NelderMeadImproveNonlinear(std::vector<cv::Point2d>& pointList, VECBITMAP<float>& dsi, VECBITMAP<float>& disp, VECBITMAP<Plane>& coeffs)
{

	int NelderMeadOptimize(float *x, int dims, float(*feval)(float*, int), int maxiters = 0);

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

	std::vector<int> xcoord(regionSize);
	std::vector<int> ycoord(regionSize);
	for (int i = 0; i < regionSize; i++) {
		xcoord[i] = pointList[i].x;
		ycoord[i] = pointList[i].y;
	}
	nth_element(xcoord.begin(), xcoord.begin() + regionSize / 2, xcoord.end());
	nth_element(ycoord.begin(), ycoord.begin() + regionSize / 2, ycoord.end());
	int x0 = *(xcoord.begin() + regionSize / 2);
	int y0 = *(ycoord.begin() + regionSize / 2);

	// Do centralization for the ease of optimization
	// If not centralized, to obtained a fairly symmetric quadratic surface my require large pertubation
	// of the value of D, E, F, which depends on the scale of current coordinates. seems not stable.
	std::vector<cv::Point2d> centralizedPointList(regionSize);
	for (int i = 0; i < regionSize; i++) {
		centralizedPointList[i].x = pointList[i].x - x0;
		centralizedPointList[i].y = pointList[i].y - y0;
	}

	//printf("finish centalization.\n");
	float a, b, c, A, B, C, D, E, F;
	a = coeffs[y0][x0].a;
	b = coeffs[y0][x0].b;
	c = coeffs[y0][x0].c;
	A = B = C = 0;
	D = a; E = b;
	F = c + D*x0 + E*y0;

	QuadraticSurface initplane(A, B, C, D, E, F);

	float vertices[6 * 7];
	int max_iters = std::min(regionSize, 15);

	//initplane.SetToArray(&vertices[0]);
	//printf("before iter.\n");
	for (int iter = 0; iter < max_iters; iter++) {

		// initialize simplex vertices
		for (int i = 1; i < 7; i++) {
			initplane.RandomPerturbTo(&vertices[6 * i]);
		}
		//for (int i = 0; i < 7; i++) {
		//	initplane.SetToArray(&vertices[6 * i]);
		//}
		// always include the initial plane 
		initplane.SetToArray(&vertices[6]);

		
		// invoke nelder-mead
		nm_opt_struct.dsi = &dsi;
		nm_opt_struct.pointList = &centralizedPointList;
		nm_opt_struct.x0 = x0;
		nm_opt_struct.y0 = y0;
		//printf("fin\n");
		float cost_before = nm_compute_quadratic_surface_cost(vertices);
		//printf("invoking..\n");
		NelderMeadOptimize(vertices, 6, nm_compute_quadratic_surface_cost, 30);
		//printf("after invoking..\n");
		float cost_after = nm_compute_quadratic_surface_cost(vertices);

		if (cost_after - cost_before > 0) {
			printf("BUG: energy increased!\n");
		}
	}

	//printf("surface neldermead done.\n");
	QuadraticSurface bestsurface(vertices[0], vertices[1], vertices[2], vertices[3], vertices[4], vertices[5]);

	for (int i = 0; i < regionSize; i++) {
		int y = centralizedPointList[i].y;
		int x = centralizedPointList[i].x;
		if (!InBound(y + y0, x + x0)) {
			printf("OUT OF BOUND BUG!\n");
		}
		disp[y + y0][x + x0] = bestsurface.ToDisparity(y, x);
	}

}

void NelderMeadEstimate(std::vector<cv::Point2d>& pointList, VECBITMAP<float>& dsi, VECBITMAP<float>& disp, VECBITMAP<Plane>& coeffs)
{

	int NelderMeadOptimize(float *x, int dims, float(*feval)(float*, int), int maxiters = 0);

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
		NelderMeadOptimize(vertices, 3, nm_compute_plane_cost, 30);
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
	//meanShiftSegmentation(imL, 2, 2, 100, g_segments);
	meanShiftSegmentation(imL, 4, 5, 20, g_segments);
	//slicSegmentation(imL, 300, 20, g_segments);
	Timer::toc();

	//Timer::tic("Intersect meanshift and SLIC");
	//cv::Mat segment1, segment2;
	//meanShiftSegmentation(imL, 4, 5, 20, segment1);
	//slicSegmentation(imL, 300, 20, segment2);
	//segment1.convertTo(segment1, CV_32FC3);
	//segment2.convertTo(segment2, CV_32FC3);
	//g_segments = (segment1 + segment2) / 2;
	//g_segments.convertTo(g_segments, CV_8UC3);
	//Timer::toc();

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

	for (int id = 0; id < nlables; id++) {
		// Optimize nonlinear part
		NelderMeadImproveNonlinear(regionList[id], dsiL, dispL, coeffsL);
	}

	
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
	//PlanefitView(imR, dsiR, coeffsR, dispR);

	//VECBITMAP<float> weightsL = PrecomputeWeights(imL);
	//VECBITMAP<float> weightsR = PrecomputeWeights(imR);

	//Timer::tic("Planefit Postprocess");
	//PostProcess(weightsL, weightsR, coeffsL, coeffsR, dispL, dispR);


	EvaluateDisparity(dispL, 0.5f, coeffsL);
}



