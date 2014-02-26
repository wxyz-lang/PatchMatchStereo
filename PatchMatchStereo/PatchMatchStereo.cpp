//
//#include <cstdio>
//#include <cstdlib>
//#include <cmath>
//#include <cstring>
//#include <ctime>
//#include <algorithm>
//#include <vector>
//#include <list>
//
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//
//#include <omp.h>
//#include "SLIC.h"
//
//#ifdef _DEBUG
//#pragma comment(lib, "opencv_core248d.lib")
//#pragma comment(lib, "opencv_highgui248d.lib")
//#pragma comment(lib, "opencv_imgproc248d.lib")
//#pragma comment(lib, "opencv_features2d248d.lib")
//#pragma comment(lib, "opencv_calib3d248d.lib")
//#pragma comment(lib, "opencv_video248d.lib")
//#pragma comment(lib, "opencv_flann248d.lib")
//#else
//#pragma comment(lib, "opencv_core248.lib")
//#pragma comment(lib, "opencv_highgui248.lib")
//#pragma comment(lib, "opencv_imgproc248.lib")
//#pragma comment(lib, "opencv_features2d248.lib")
//#pragma comment(lib, "opencv_calib3d248.lib")
//#pragma comment(lib, "opencv_video248.lib")
//#pragma comment(lib, "opencv_flann248.lib")
//#endif
//
//
//#define USE_OPENMP
////#define FRONTAL_PARALLEL_ONLY
////#define LOAD_RESULT_FROM_LAST_RUN
//#define DO_POST_PROCESSING
////#define USE_SUBPIXEL_MATCHING_COST
////#define USE_NELDERMEAD_OPT
//
//
//
//// Gobal variables
//static int nrows, ncols;
//static float BAD_PLANE_PENALTY = 120;  // defined as 2 times the max cost of dsi.
//static const float	gamma_proximity = 25;
//static int			g_improve_cnt = 0;
//
//static const int	patch_w		= 35;
//static const int	patch_r		= 17;
//static const int	maxiters	= 2;
//static const float  alpha		= 0.9;
//static const float	gamma		= 10;
//static const double	tau_col		= 10;
//static const double tau_grad	= 2;
//
//static const int folder_id = 2;
//static const std::string folders[] = { "tsukuba/", "venus/", "teddy/", "cones/", "Bowling2/", "Baby1/" };
//static const int scales[]	= { 16, 8, 4, 4, 3, 3 };
//static const int drange[]	= { 16, 20, 60, 60, 70, 70 };
//static const int scale		= scales[folder_id];
//static const int ndisps		= drange[folder_id];
//static const int dmax		= ndisps - 1;
//
//cv::Mat g_L, g_OC, g_ALL, g_DISC, g_GT, g_segments;
//cv::Mat g_stereo, g_mydisp;
//
//
//using namespace std;
//using namespace cv;
//
//
//struct Plane {
//	float a, b, c;
//	float nx, ny, nz;
//	Plane() {}
//	Plane(float a_, float b_, float c_, float nx_, float ny_, float nz_)
//		:a(a_), b(b_), c(c_), nx(nx_), ny(ny_), nz(nz_) {}
//	Plane(float nx_, float ny_, float nz_, int y, int x, float z)
//	{
//		const float eps = 0.001;
//		nx = nx_;  ny = ny_;  nz = nz_;
//		if (std::abs(nz) < eps) {
//			if (nz > 0)  nz = +eps;
//			else		 nz = -eps;
//		}
//		a = -nx / nz;
//		b = -ny / nz;
//		c = (nx * x + ny * y + nz * z) / nz;
//	}
//	Plane ReparametrizeInOtherView(int y, int x, int sign, int& qy, int &qx)
//	{
//		float z = a * x + b * y + c;
//		qx = x + sign * z;
//		qy = y;
//		return Plane(nx, ny, nz, qy, qx, z);
//	}
//	void RandomAssign()
//	{
//		const int RAND_HALF = RAND_MAX / 2;
//		float y = nrows / 2.f, x = ncols / 2.f;
//		float z = dmax * ((double)rand() / RAND_MAX);
//		float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
//		float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
//		float nz = ((double)rand() - RAND_HALF) / RAND_HALF;
//		float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
//		nx /= norm;
//		ny /= norm;
//		nz /= norm;
//		*this = Plane(nx, ny, nz, y, x, z);
//	}
//	Plane RandomSearch(int y, int x, float radius_z0, float radius_n)
//	{
//		const int RAND_HALF = RAND_MAX / 2;
//
//		float nx__ = nx + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);
//		float ny__ = ny + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);
//		float nz__ = nz + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);
//
//		float z = a * x + b * y + c
//			+ radius_z0 * (((double)rand() - RAND_HALF) / RAND_HALF);
//		z = std::max(0.f, std::min(z, (float)dmax));
//		float norm = std::max(0.01f, sqrt(nx__*nx__ + ny__*ny__ + nz__*nz__));
//		nx__ /= norm;
//		ny__ /= norm;
//		nz__ /= norm;
//#ifndef FRONTAL_PARALLEL_ONLY
//		return Plane(nx__, ny__, nz__, y, x, z);
//#else
//		z = (int)(z + 0.5);
//		return Plane(0.f, 0.f, 1.f, y, x, z);
//#endif
//	}
//	float ToDisparity(int y, int x) { return a * x + b * y + c; }
//	void GetAbc(float *v) { v[0] = a; v[1] = b; v[2] = c; }
//	void SetAbc(float *abc)
//	{
//		a = abc[0];
//		b = abc[1];
//		c = abc[2];
//	}
//};
//
//template<class T>
//class VECBITMAP {
//public:
//	T *data;
//	int w, h, n;
//	bool is_shared;
//	T *get(int y, int x) { return &data[(y*w + x)*n]; }		/* Get patch (y, x). */
//	T *line_n1(int y) { return &data[y*w]; }				/* Get line y assuming n=1. */
//	VECBITMAP() { data = NULL; }
//	VECBITMAP(const VECBITMAP& obj)
//	{
//		// This constructor is very necessary for returning an object in a function,
//		// in the case that the Name Return Value Optimization (NRVO) is turned off.
//		w = obj.w; h = obj.h; n = obj.n; is_shared = obj.is_shared;
//		if (is_shared) { data = obj.data; }
//		else { data = new T[w*h*n]; memcpy(data, obj.data, w*h*n*sizeof(T)); }
//
//	}
//	VECBITMAP(int h_, int w_, int n_ = 1, T* data_ = NULL)
//	{
//		w = w_; h = h_; n = n_;
//		if (!data_) { data = new T[w*h*n]; is_shared = false; }
//		else	    { data = data_;        is_shared = true; }
//	}
//	VECBITMAP& operator=(const VECBITMAP& m)
//	{
//		if (data) { delete[] data; }
//		w = m.w; h = m.h; n = m.n; is_shared = m.is_shared;
//		if (m.is_shared) { data = m.data; }
//		else { data = new T[w*h*n]; memcpy(data, m.data, w*h*n*sizeof(T)); }
//		return *this;
//	}
//	~VECBITMAP() { if (!is_shared) delete[] data; }
//	T *operator[](int y) { return &data[y*w]; }
//
//	void SaveToBinaryFile(std::string filename)
//	{
//		FILE *fid = fopen(filename.c_str(), "wb");
//		assert(fid != NULL);
//		fwrite(data, sizeof(T), w*h*n, fid);
//		fclose(fid);
//	}
//	void LoadFromBinaryFile(std::string filename)
//	{
//		FILE *fid = fopen(filename.c_str(), "rb");
//		assert(fid != NULL);
//		fread(data, sizeof(T), w*h*n, fid);
//		fclose(fid);
//	}
//};
//
//VECBITMAP<Plane> g_coeffsL, g_coeffsR;
//VECBITMAP<float> g_colgradL, g_colgradR;
//VECBITMAP<unsigned char> g_imL, g_imR;
//
//
//inline bool InBound(float y, float x) { return 0 <= y && y < nrows && 0 <= x && x < ncols; }
//
//VECBITMAP<float> ComputeLocalPatchWeights(int yc, int xc, VECBITMAP<unsigned char>& img)
//{
//	VECBITMAP<float> ret(patch_w, patch_w);
//	float *w = ret.data;
//	memset(w, 0, patch_w * patch_w * sizeof(float));
//	unsigned char *rgb1 = img.get(yc, xc);
//	unsigned char *rgb2 = NULL;
//
//	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
//	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);
//
//	#pragma omp parallel for
//	for (int y = yb; y <= ye; y++) {
//		for (int x = xb; x <= xe; x++) {
//			rgb2 = img.get(y, x);
//			float cost = std::abs((float)rgb1[0] - (float)rgb2[0])
//				+ std::abs((float)rgb1[1] - (float)rgb2[1])
//				+ std::abs((float)rgb1[2] - (float)rgb2[2]);
//			w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] = cost;
//		}
//	}
//	for (int i = 0; i < patch_w * patch_w; i++) {
//		w[i] = exp(-w[i] / gamma);
//	}
//	return ret;
//}
//
//float *PrecomputeWeights(VECBITMAP<unsigned char>& im)
//{
//	const int nw = nrows * ncols * patch_w * patch_w;
//	const int patchsize = patch_w * patch_w;
//	float *weights = new float[nw];
//	memset(weights, 0, nw * sizeof(float));
//
//	float w_prox[patch_w * patch_w];
//	for (int y = 0; y < patch_w; y++) {
//		for (int x = 0; x < patch_w; x++) {
//			float dist = sqrt((y - patch_r) * (y - patch_r) + (x - patch_r) * (x - patch_r));
//			w_prox[y * patch_w + x] = exp(-dist / gamma_proximity);
//		}
//	}
//
//	#pragma omp parallel for
//	for (int yc = 0; yc < nrows; yc++) {
//		for (int xc = 0; xc < ncols; xc++) {
//
//			float *w = &weights[(yc*ncols + xc) * patchsize];
//			unsigned char *rgb1 = im.get(yc, xc);
//			unsigned char *rgb2 = NULL;
//
//			int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
//			int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);
//
//			for (int y = yb; y <= ye; y++) {
//				for (int x = xb; x <= xe; x++) {
//					rgb2 = im.get(y, x);
//					float cost = std::abs((float)rgb1[0] - (float)rgb2[0])
//						+ std::abs((float)rgb1[1] - (float)rgb2[1])
//						+ std::abs((float)rgb1[2] - (float)rgb2[2]);
//					w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] = cost;
//				}
//			}
//		}
//	}
//
//	for (int i = 0; i < nw; i++) {
//		weights[i] = exp(-weights[i] / gamma);
//	}
//
//	return weights;
//}
//
//double ComputePlaneCost(int yc, int xc, Plane& coeff_try, VECBITMAP<float>& dsi, VECBITMAP<float>& w)
//{
//	double cost = 0;
//	for (int y = yc - patch_r; y <= yc + patch_r; y++) {
//		for (int x = xc - patch_r; x <= xc + patch_r; x++) {
//			int d = 0.5 + (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
//			if (InBound(y, x)) {
//				if (d < 0 || d > dmax) {	// must be a bad plane.
//					cost += BAD_PLANE_PENALTY;
//				}
//				else {
//					cost += w[y - yc + patch_r][x - xc + patch_r] * (dsi.get(y, x)[d]);
//				}
//			}
//		}
//	}
//	return cost;
//}
//
//double ComputePlaneCost(int yc, int xc, Plane& coeff_try, VECBITMAP<float>& imL, VECBITMAP<float> &imR, VECBITMAP<float>& w, int sign)
//{
//	BAD_PLANE_PENALTY = 3 * ((1 - alpha) * tau_col + alpha * tau_grad);		// set to 3 times the max possible matchin cost.
//	double cost = 0;
//	for (int y = yc - patch_r; y <= yc + patch_r; y++) {
//		for (int x = xc - patch_r; x <= xc + patch_r; x++) {
//			//int d = 0.5 + (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
//			float d = (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
//			if (InBound(y, x)) {
//				
//
//				if (d < 0 || d > dmax) {	// must be a bad plane.
//					cost += BAD_PLANE_PENALTY;
//				}
//				else {
//					//int xm = (x + sign*d + ncols) % ncols;
//					//float *pL = imL.get(y, x);
//					//float *pR = imR.get(y, xm);
//					//float cost_col = fabs(pL[0] - pR[0]) + fabs(pL[1] - pR[1]) + fabs(pL[2] - pR[2]);
//					//float cost_grad = fabs(pL[3] - pR[3]) + fabs(pL[4] - pR[4]);
//					//cost_col = std::min((float)tau_col, cost_col);
//					//cost_grad = std::min((float)tau_grad, cost_grad);
//					//double subcost = (1 - alpha) * cost_col + alpha * cost_grad;
//					//cost += w[y - yc + patch_r][x - xc + patch_r] * subcost;
//
//
//					float xm, xmL, xmR, wL, wR;
//					double costL, costR, cost_col, cost_grad;
//					xm	= std::max(0.f, std::min((float)ncols - 1, x + sign * d));
//					xmL = (int)(xm);
//					xmR = (int)(xm + 0.99);
//					wL	= xmR - xm;
//					wR	= 1 - wL;
//
//					float *pL = imL.get(y, x);
//					float *pRmL = imR.get(y, xmL);
//					float *pRmR = imR.get(y, xmR);
//					float pR[5];
//					for (int i = 0; i < 5; i++) {
//						pR[i] = wL * pRmL[i] + wR * pRmR[i];
//					}
//					//float *pR = imR.get(y, (int)(xm + 0.5));
//					cost_col = fabs(pL[0] - pR[0])
//							 + fabs(pL[1] - pR[1])
//							 + fabs(pL[2] - pR[2]);
//					cost_col = std::min(tau_col, cost_col);
//					cost_grad = fabs(pL[3] - pR[3])
//							  + fabs(pL[4] - pR[4]);
//					cost_grad = std::min(tau_grad, cost_grad);
//					cost += w[y - yc + patch_r][x - xc + patch_r] * ((1 - alpha) * cost_col + alpha * cost_grad);
//				}
//			}
//		}
//	}
//	return cost;
//}
//
//void RandomInit(VECBITMAP<Plane>& coeffs, VECBITMAP<double>& bestcosts, VECBITMAP<float>& dsi, float* weights)
//{
//	const int	RAND_HALF = RAND_MAX / 2;
//	const float eps = 0.001;
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			// FIXME: Consider using double to avoid numerical issues
//			
//#ifndef FRONTAL_PARALLEL_ONLY
//			float z = dmax * ((double)rand() / RAND_MAX);
//			float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float nz = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
//			nx /= norm;
//			ny /= norm;
//			nz /= norm;
//#else 
//			//float z = (int)((4 * ndisps - 1) * ((double)rand() / RAND_MAX) + 0.5);
//			//z /= 4.f;
//			float z = (int)(dmax * ((double)rand() / RAND_MAX));
//			float nx = 0;
//			float ny = 0;
//			float nz = 1.f;
//#endif
//
//			coeffs[y][x] = Plane(nx, ny, nz, y, x, z);
//		}
//	}
//
//	#pragma omp parallel for
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			bestcosts[y][x] = ComputePlaneCost(y, x, coeffs[y][x], dsi, VECBITMAP<float>(patch_w, patch_w, 1, &weights[(y*ncols + x) * patch_w * patch_w]));
//		}
//	}
//}
//
//void RandomInit(VECBITMAP<Plane>& coeffs, VECBITMAP<double>& bestcosts, VECBITMAP<float>& colgradL, VECBITMAP<float>& colgradR, float* weights, int sign)
//{
//	const int	RAND_HALF = RAND_MAX / 2;
//	const float eps = 0.001;
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			// FIXME: Consider using double to avoid numerical issues
//
//#ifndef FRONTAL_PARALLEL_ONLY
//			float z = dmax * ((double)rand() / RAND_MAX);
//			float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float nz = ((double)rand() - RAND_HALF) / RAND_HALF;
//			float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
//			nx /= norm;
//			ny /= norm;
//			nz /= norm;
//#else 
//			//float z = (int)((4 * ndisps - 1) * ((double)rand() / RAND_MAX) + 0.5);
//			//z /= 4.f;
//			float z = (int)(dmax * ((double)rand() / RAND_MAX));
//			float nx = 0;
//			float ny = 0;
//			float nz = 1.f;
//#endif
//
//			coeffs[y][x] = Plane(nx, ny, nz, y, x, z);
//		}
//	}
//
//	#pragma omp parallel for
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			bestcosts[y][x] = ComputePlaneCost(y, x, coeffs[y][x], colgradL, colgradR, 
//				VECBITMAP<float>(patch_w, patch_w, 1, &weights[(y*ncols + x) * patch_w * patch_w]), sign);
//		}
//	}
//}
//
//void RandomAbc(float *v, int y, int x)
//{
//	const int	RAND_HALF = RAND_MAX / 2;
//	const float eps = 0.001;
//
//#ifndef FRONTAL_PARALLEL_ONLY
//	float z = dmax * ((double)rand() / RAND_MAX);
//	float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
//	float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
//	float nz = ((double)rand() - RAND_HALF) / RAND_HALF;
//	float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
//	nx /= norm;
//	ny /= norm;
//	nz /= norm;
//#else
//	//float z = (int)((4 * ndisps - 1) * ((double)rand() / RAND_MAX) + 0.5);
//	//z /= 4.f;
//	float z = (int)(dmax * ((double)rand() / RAND_MAX));
//	float nx = 0;
//	float ny = 0;
//	float nz = 1.f;
//#endif
//
//	Plane(nx, ny, nz, y, x, z).GetAbc(v);
//}
//
//void ImproveGuess(int y, int x, Plane& coeff_old, double& bestcost, Plane& coeff_try, VECBITMAP<float>& dsi, VECBITMAP<float>& w)
//{
//	float cost = ComputePlaneCost(y, x, coeff_try, dsi, w);
//	if (cost < bestcost) {
//		g_improve_cnt++;
//		bestcost = cost;
//		coeff_old = coeff_try;
//	}
//}
//
//void ImproveGuess(int y, int x, Plane& coeff_old, double& bestcost, Plane& coeff_try, VECBITMAP<float>& imL, VECBITMAP<float>& imR, VECBITMAP<float>& w, int sign)
//{
//	float cost = ComputePlaneCost(y, x, coeff_try, imL, imR, w, sign);
//	if (cost < bestcost) {
//		g_improve_cnt++;
//		bestcost = cost;
//		coeff_old = coeff_try;
//	}
//}
//
//struct NM_OPT_PARAM {
//	int yc, xc;
//	Plane coeff;
//	VECBITMAP<float> *dsi;
//	VECBITMAP<float> *w;
//	VECBITMAP<float> *colgradL;
//	VECBITMAP<float> *colgradR;
//	int sign;
//};
//NM_OPT_PARAM nm_opt_struct;
//
//int NelderMeadOptimize(float *, float(*)(float*, int), int);
//
//float nm_compute_plane_cost(float *abc, int n = 3)
//{
//	Plane coeff;
//	coeff.SetAbc(abc);
//	return ComputePlaneCost(nm_opt_struct.yc, nm_opt_struct.xc, coeff, *nm_opt_struct.dsi, *nm_opt_struct.w);
//}
//
//float nm_compute_plane_cost_subpixel(float *abc, int n = 3)
//{
//	Plane coeff;
//	coeff.SetAbc(abc);
//	return ComputePlaneCost(nm_opt_struct.yc, nm_opt_struct.xc, coeff, *nm_opt_struct.colgradL, *nm_opt_struct.colgradR, *nm_opt_struct.w, nm_opt_struct.sign);
//}
//
//void ProcessViewDelegate(int y, int x,
//	VECBITMAP<Plane>& coeffsL,		VECBITMAP<Plane>& coeffsR,
//	VECBITMAP<double>& bestcostsL,	VECBITMAP<double>& bestcostsR,
//	VECBITMAP<float>& dsiL,			VECBITMAP<float>& dsiR,
//	float* weightsL,				float* weightsR,
//	int iter,
//	int sign)
//{
//	int xchange, ychange;
//	if (iter % 2 == 0)  xchange = ychange = +1;
//	else				xchange = ychange = -1;
//	const int dx[] = { -1, 0, +1, 0 };
//	const int dy[] = { 0, -1, 0, +1 };
//	float nm_opt_x[3 * 4];
//
//
//	VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y*ncols + x) * patch_w * patch_w]);
//#ifndef USE_NELDERMEAD_OPT
//	// Spatial Propagation
//	int qy = y - ychange, qx = x;
//	if (InBound(qy, qx)) {
//		Plane coeff_try = coeffsL[qy][qx];
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
//	}
//
//	qy = y; qx = x - xchange;
//	if (InBound(qy, qx)) {
//		Plane coeff_try = coeffsL[qy][qx];
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
//	}
//	//int qx, qy;
//	//for (int dir = 0; dir < 4; dir++) {
//	//	qy = y + dy[dir];
//	//	qx = x + dx[dir];
//	//	if (InBound(qy, qx)) {
//	//		Plane coeff_try = coeffsL[qy][qx];
//	//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
//	//	}
//	//}
//
//	// Random Search
//	float radius_z = dmax / 2.0f;
//	float radius_n = 1.0f;
//	while (radius_z >= 0.1) {
//		Plane coeff_try = coeffsL[y][x].RandomSearch(y, x, radius_z, radius_n);
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
//		radius_z /= 2.0f;
//		radius_n /= 2.0f;
//	}
//
//	// View Propagation
//	Plane coeff_try = coeffsL[y][x].ReparametrizeInOtherView(y, x, sign, qy, qx);
//	if (0 <= qx && qx < ncols) {
//		VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(qy*ncols + qx) * patch_w * patch_w]);
//		ImproveGuess(qy, qx, coeffsR[qy][qx], bestcostsR[qy][qx], coeff_try, dsiR, wR);
//	}
//#else
//
//	for (int dir = 0; dir < 4; dir++) {
//		int qy = y + dy[dir];
//		int qx = x + dx[dir];
//		if (InBound(qy, qx)) {
//			coeffsL[qy][qx].GetAbc(nm_opt_x + 3 * dir);
//		}
//		else {
//			RandomAbc(nm_opt_x + 3 * dir, qy, qx);
//		}
//	}
//	nm_opt_struct.yc = y;
//	nm_opt_struct.xc = x;
//	nm_opt_struct.w = &wL;
//	nm_opt_struct.dsi = &dsiL;
//	NelderMeadOptimize(nm_opt_x, nm_compute_plane_cost, 5);
//	coeffsL[y][x].SetAbc(nm_opt_x);
//#endif
//}
//
//void ProcessView(
//	VECBITMAP<Plane>& coeffsL,		VECBITMAP<Plane>& coeffsR, 
//	VECBITMAP<double>& bestcostsL,	VECBITMAP<double>& bestcostsR,
//	VECBITMAP<float>& dsiL,			VECBITMAP<float>& dsiR, 
//	float* weightsL,				float* weightsR,
//	int iter, 
//	int sign)
//{
//	if (iter % 2 == 0) {
//		#pragma omp parallel for
//		for (int y = 0; y < nrows; y++) {
//			for (int x = 0; x < ncols; x++) {
//				ProcessViewDelegate(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, dsiL, dsiR, weightsL, weightsR, iter, sign);
//			}
//		}
//	}
//	else {
//		#pragma omp parallel for
//		for (int y = nrows - 1; y >= 0; y--) {
//			for (int x = ncols - 1; x >= 0; x--) {
//				ProcessViewDelegate(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, dsiL, dsiR, weightsL, weightsR, iter, sign);
//			}
//		}
//	}
//}
//
//void ProcessViewSubpixelDelegate(int y, int x,
//	VECBITMAP<Plane>& coeffsL,		VECBITMAP<Plane>& coeffsR,
//	VECBITMAP<double>& bestcostsL,	VECBITMAP<double>& bestcostsR,
//	VECBITMAP<float>& colgradL,		VECBITMAP<float>& colgradR,
//	float* weightsL,				float* weightsR,
//	int iter,
//	int sign)
//{
//	int xchange, ychange;
//	if (iter % 2 == 0)  xchange = ychange = +1;
//	else				xchange = ychange = -1;
//
//	const int dx[] = { -1, 0, +1, 0 };
//	const int dy[] = { 0, -1, 0, +1 };
//	float nm_opt_x[3 * 4];
//
//
//	VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y*ncols + x) * patch_w * patch_w]);
//#ifndef USE_NELDERMEAD_OPT
//	// Spatial Propagation
//	int qy = y - ychange, qx = x;
//	if (InBound(qy, qx)) {
//		Plane coeff_try = coeffsL[qy][qx];
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, colgradL, colgradR, wL, sign);
//	}
//
//	qy = y; qx = x - xchange;
//	if (InBound(qy, qx)) {
//		Plane coeff_try = coeffsL[qy][qx];
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, colgradL, colgradR, wL, sign);
//	}
//
//	//int qx, qy;
//	//for (int dir = 0; dir < 4; dir++) {
//	//	qy = y + dy[dir];
//	//	qx = x + dx[dir];
//	//	if (InBound(qy, qx)) {
//	//		Plane coeff_try = coeffsL[qy][qx];
//	//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, colgradL, colgradR, wL, sign);
//	//	}
//	//}
//
//	// Random Search
//	float radius_z = dmax / 2.0f;
//	float radius_n = 1.0f;
//	while (radius_z >= 0.1) {
//		Plane coeff_try = coeffsL[y][x].RandomSearch(y, x, radius_z, radius_n);
//		ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, colgradL, colgradR, wL, sign);
//		radius_z /= 2.0f;
//		radius_n /= 2.0f;
//	}
//
//	// View Propagation
//	Plane coeff_try = coeffsL[y][x].ReparametrizeInOtherView(y, x, sign, qy, qx);
//	if (0 <= qx && qx < ncols) {
//		VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(qy*ncols + qx) * patch_w * patch_w]);
//		ImproveGuess(qy, qx, coeffsR[qy][qx], bestcostsR[qy][qx], coeff_try, colgradR, colgradL, wR, -sign);
//	}
//#else
//
//	for (int dir = 0; dir < 4; dir++) {
//		int qy = y + dy[dir];
//		int qx = x + dx[dir];
//		if (InBound(qy, qx)) {	// FIXME: only propagate from UP and left, RANDOM initialize other two vertices.
//			coeffsL[qy][qx].GetAbc(nm_opt_x + 3 * dir);
//		}
//		else {
//			RandomAbc(nm_opt_x + 3 * dir, qy, qx);
//		}
//	}
//
//	nm_opt_struct.yc = y;
//	nm_opt_struct.xc = x;
//	nm_opt_struct.w = &wL;
//	nm_opt_struct.sign = sign;
//
//	if (sign == -1) {
//		nm_opt_struct.colgradL = &colgradL;
//		nm_opt_struct.colgradR = &colgradR;
//	}
//	else {
//		nm_opt_struct.colgradL = &colgradR;
//		nm_opt_struct.colgradR = &colgradL;
//	}
//
//	NelderMeadOptimize(nm_opt_x, nm_compute_plane_cost_subpixel, 5);
//	coeffsL[y][x].SetAbc(nm_opt_x);
//#endif
//}
//
//void ProcessViewSubpixel(
//	VECBITMAP<Plane>& coeffsL,		VECBITMAP<Plane>& coeffsR,
//	VECBITMAP<double>& bestcostsL,	VECBITMAP<double>& bestcostsR,
//	VECBITMAP<float>& colgradL,		VECBITMAP<float>& colgradR,
//	float* weightsL, float* weightsR,
//	int iter,
//	int sign)
//{
//	if (iter % 2 == 0) {
//		#pragma omp parallel for
//		for (int y = 0; y < nrows; y++) {
//			for (int x = 0; x < ncols; x++) {
//				ProcessViewSubpixelDelegate(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, colgradL, colgradR, weightsL, weightsR, iter, sign);
//			}
//		}
//	}
//	else {
//		#pragma omp parallel for
//		for (int y = nrows - 1; y >= 0; y--) {
//			for (int x = ncols - 1; x >= 0; x--) {
//				ProcessViewSubpixelDelegate(y, x, coeffsL, coeffsR, bestcostsL, bestcostsR, colgradL, colgradR, weightsL, weightsR, iter, sign);
//			}
//		}
//	}
//}
//
//void WinnerTakesAll(VECBITMAP<float>& dsi, VECBITMAP<float> &disp)
//{
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			int minidx = 0;
//			float *cost = dsi.get(y, x);
//			for (int k = 1; k < ndisps; k++) {
//				if (cost[k] < cost[minidx]) {
//					minidx = k;
//				}
//			}
//			disp[y][x] = minidx;
//		}
//	}
//}
//
//void PlaneMapToDisparityMap(VECBITMAP<Plane>& coeffs, VECBITMAP<float>& disp)
//{
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			disp[y][x] = coeffs[y][x].ToDisparity(y, x);
//		}
//	}
//}
//
//void CrossCheck(VECBITMAP<float>& dispL, VECBITMAP<float>& dispR, VECBITMAP<bool>& validL, VECBITMAP<bool>& validR)
//{
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//
//			int xR = std::max(0.f, std::min((float)ncols, x - dispL[y][x]));
//			validL[y][x] = (std::abs(dispL[y][x] - dispR[y][xR]) <= 1);
//
//			int xL = std::max(0.f, std::min((float)ncols, x + dispR[y][x]));
//			validR[y][x] = (std::abs(dispR[y][x] - dispL[y][xL]) <= 1);
//		}
//	}
//}
//
//void FillHole(int y, int x, VECBITMAP<bool>& valid, VECBITMAP<Plane>& coeffs)
//{
//	// This function fills the invalid pixel (y,x) by finding its nearst (left and right) 
//    // valid neighbors on the same scanline, and select the one with lower disparity.
//
//	int xL = x - 1, xR = x + 1, bestx = x;
//	while (!valid[y][xL] && 0 <= xL) {
//		xL--;
//	}
//	while (!valid[y][xR] && xR < ncols) {
//		xR++;
//	}
//	if (0 <= xL) {
//		bestx = xL;
//	}
//	if (xR < ncols) {
//		if (bestx == xL) {
//			float dL = coeffs[y][xL].ToDisparity(y, x);
//			float dR = coeffs[y][xR].ToDisparity(y, x);
//			if (dR < dL) {
//				bestx = xR;
//			}
//		}
//		else {
//			bestx = xR;
//		}
//	}
//	coeffs[y][x] = coeffs[y][bestx];
//}
//
//void WeightedMedianFilter(int yc, int xc, VECBITMAP<float>& disp, VECBITMAP<float>& weights, VECBITMAP<bool>& valid, bool useInvalidPixels)
//{
//	std::vector<std::pair<float, float>> dw_pairs;
//
//	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
//	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);
//
//	for (int y = yb; y <= ye; y++) {
//		for (int x = xb; x <= xe; x++) {
//			if (useInvalidPixels || valid[y][x]) {
//				std::pair<float, float> dw(disp[y][x], weights[y - yc + patch_r][x - xc + patch_r]);
//				dw_pairs.push_back(dw);
//			}
//		}
//	}
//
//	std::sort(dw_pairs.begin(), dw_pairs.end());
//
//	float w = 0.f, wsum = 0.f;
//	for (int i = 0; i < dw_pairs.size(); i++) {
//		wsum += dw_pairs[i].second;
//	}
//
//	for (int i = 0; i < dw_pairs.size(); i++) {
//		w += dw_pairs[i].second;
//		if (w >= wsum / 2.f) {
//			// Note that this line can always be reached.
//			if (i > 0) {
//				disp[yc][xc] = (dw_pairs[i - 1].first + dw_pairs[i].first) / 2.f;
//			}
//			else {
//				disp[yc][xc] = dw_pairs[i].first;
//			}
//			break;
//		}
//	}
//}
//
//void PostProcess(
//	float *	weightsL,			float *weightsR, 
//	VECBITMAP<Plane>& coeffsL,	VECBITMAP<Plane>& coeffsR, 
//	VECBITMAP<float>& dispL,	VECBITMAP<float>& dispR)
//{
//	// This function perform several times of weighted median filtering at each invalid position.
//	// weights of the neighborhood are set by exp(-||cp-cq|| / gamma), except in the last iteration,
//	// where the weights of invalid pixels are set to zero.
//
//	VECBITMAP<bool> validL(nrows, ncols), validR(nrows, ncols);
//	PlaneMapToDisparityMap(coeffsL, dispL);
//	PlaneMapToDisparityMap(coeffsR, dispR);
//#ifdef DO_POST_PROCESSING
//	// Hole filling
//	CrossCheck(dispL, dispR, validL, validR);
//	#pragma omp parallel for
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			if (!validL[y][x]) {
//				FillHole(y, x, validL, coeffsL);
//			}
//			if (!validR[y][x]) {
//				FillHole(y, x, validR, coeffsR);
//			}
//		}
//	}
//
//	// Weighted median filtering 
//	int maxround = 1; 
//	bool useInvalidPixels = true;
//	for (int round = 0; round < maxround; round++) {
//
//		PlaneMapToDisparityMap(coeffsL, dispL);
//		PlaneMapToDisparityMap(coeffsR, dispR);
//		CrossCheck(dispL, dispR, validL, validR);
//
//		//if (round + 1 == maxround) {
//		//	useInvalidPixels = false;
//		//}
//		#pragma omp parallel for
//		for (int y = 0; y < nrows; y++) {
//			if (y % 10 == 0) { printf("median filtering row %d\n", y); }
//			for (int x = 0; x < ncols; x++) {
//				if (!validL[y][x]){
//					VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y * ncols + x) * patch_w * patch_w]);
//					WeightedMedianFilter(y, x, dispL, wL, validL, useInvalidPixels);
//				}
//				if (!validR[y][x]){
//					VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(y * ncols + x) * patch_w * patch_w]);
//					WeightedMedianFilter(y, x, dispR, wR, validR, useInvalidPixels);
//				}
//			}
//		}
//	}
//#endif
//}
//
//void RunPatchMatchStereo(VECBITMAP<unsigned char>& imL, VECBITMAP<unsigned char>& imR, 
//						 VECBITMAP<float>& dsiL,		VECBITMAP<float>& dsiR,  
//						 VECBITMAP<float>& dispL,		VECBITMAP<float>& dispR)
//{
//	VECBITMAP<Plane> coeffsL(nrows, ncols), coeffsR(nrows, ncols);
//	VECBITMAP<double> bestcostsL(nrows, ncols), bestcostsR(nrows, ncols);
//
//
//	int tic, toc;
//	// Precompute color weights to speed up
//	tic = clock();  printf("precomputing weights ...\n");
//	float *weightsL = PrecomputeWeights(imL);
//	float *weightsR = PrecomputeWeights(imR);
//	toc = clock();  printf("precomputing weights use %.2fs\n\n", (toc - tic) / 1000.f);
//
//
//	// Random initialization
//	tic = clock();  printf("performing random init ...\n");
//	RandomInit(coeffsL, bestcostsL, dsiL, weightsL);
//	RandomInit(coeffsR, bestcostsR, dsiR, weightsR);
//	toc = clock();  printf("random init use %.2fs\n\n", (toc - tic) / 1000.f);
//
//#ifndef LOAD_RESULT_FROM_LAST_RUN
//	// Iteration
//	for (int iter = 0; iter < maxiters; iter++) {
//
//		tic = clock();  printf("iter %d, scaning left view ...\n", iter);
//		ProcessView(coeffsL, coeffsR, bestcostsL, bestcostsR, dsiL, dsiR, weightsL, weightsR, iter, -1);
//		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//
//		tic = clock();  printf("iter %d, scaning right view ...\n", iter);
//		ProcessView(coeffsR, coeffsL, bestcostsR, bestcostsL, dsiR, dsiL, weightsR, weightsL, iter, +1);
//		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//	}
//	printf("g_improve_cnt: %d\n", g_improve_cnt);
//
//	coeffsL.SaveToBinaryFile(folders[folder_id] + "coeffsL.bin");
//	coeffsR.SaveToBinaryFile(folders[folder_id] + "coeffsR.bin");
//#else
//	coeffsL.LoadFromBinaryFile(folders[folder_id] + "coeffsL.bin");
//	coeffsR.LoadFromBinaryFile(folders[folder_id] + "coeffsR.bin");
//#endif
//
//	// Post processing
//	tic = clock();  printf("POST PROCESSING...\n");
//	PostProcess(weightsL, weightsR, coeffsL, coeffsR, dispL, dispR);
//	toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//
//	g_coeffsL = coeffsL;
//	if (!weightsL) delete[] weightsL;
//	if (!weightsR) delete[] weightsR;
//}
//
//void RunPatchMatchStereoSubpixel(VECBITMAP<unsigned char>& imL, VECBITMAP<unsigned char>& imR,
//	VECBITMAP<float>& colgradL, VECBITMAP<float>& colgradR,
//	VECBITMAP<float>& dispL, VECBITMAP<float>& dispR)
//{
//	VECBITMAP<Plane> coeffsL(nrows, ncols), coeffsR(nrows, ncols);
//	VECBITMAP<double> bestcostsL(nrows, ncols), bestcostsR(nrows, ncols);
//
//
//	int tic, toc;
//	// Precompute color weights to speed up
//	tic = clock();  printf("precomputing weights ...\n");
//	float *weightsL = PrecomputeWeights(imL);
//	float *weightsR = PrecomputeWeights(imR);
//	toc = clock();  printf("precomputing weights use %.2fs\n\n", (toc - tic) / 1000.f);
//
//
//	// Random initialization
//	tic = clock();  printf("performing random init ...\n");
//	RandomInit(coeffsL, bestcostsL, colgradL, colgradR, weightsL, -1);
//	RandomInit(coeffsR, bestcostsR, colgradR, colgradL, weightsR, +1);
//	toc = clock();  printf("random init use %.2fs\n\n", (toc - tic) / 1000.f);
//
//#ifndef LOAD_RESULT_FROM_LAST_RUN
//	// Iteration
//	for (int iter = 0; iter < maxiters; iter++) {
//
//		tic = clock();  printf("iter %d, scaning left view ...\n", iter);
//		ProcessViewSubpixel(coeffsL, coeffsR, bestcostsL, bestcostsR, colgradL, colgradR, weightsL, weightsR, iter, -1);
//		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//
//		tic = clock();  printf("iter %d, scaning right view ...\n", iter);
//		ProcessViewSubpixel(coeffsR, coeffsL, bestcostsR, bestcostsL, colgradR, colgradL, weightsR, weightsL, iter, +1);
//		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//	}
//
//	printf("g_improve_cnt: %d\n", g_improve_cnt);
//
//	coeffsL.SaveToBinaryFile(folders[folder_id] + "coeffsL.bin");
//	coeffsR.SaveToBinaryFile(folders[folder_id] + "coeffsR.bin");
//#else
//	coeffsL.LoadFromBinaryFile(folders[folder_id] + "coeffsL.bin");
//	coeffsR.LoadFromBinaryFile(folders[folder_id] + "coeffsR.bin");
//#endif
//
//	// Post processing
//	tic = clock();  printf("POST PROCESSING...\n");
//	PostProcess(weightsL, weightsR, coeffsL, coeffsR, dispL, dispR);
//	toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
//
//	g_coeffsL = coeffsL;
//	if (!weightsL) delete[] weightsL;
//	if (!weightsR) delete[] weightsR;
//}
//
//unsigned hamdist(long long x, long long  y)
//{
//	int dist = 0;
//	long long val = x ^ y; 
//	while (val) {
//		++dist;
//		val &= val - 1;
//	}
//	return dist;
//}
//
//VECBITMAP<long long> ComputeCensusImage(VECBITMAP<unsigned char> im)
//{
//	int vpad = 3, hpad = 4;
//	VECBITMAP<long long> census(nrows, ncols);
//	memset(census.data, 0, nrows * ncols * sizeof(long long));
//
//	for (int yc = 0; yc < nrows; yc++) {
//		for (int xc = 0; xc < ncols; xc++) {
//
//			int uu = std::max(yc - vpad, 0);
//			int dd = std::min(yc + vpad, nrows - 1);
//			int ll = std::max(xc - hpad, 0);
//			int rr = std::min(xc + hpad, ncols - 1);
//
//			int idx = 0;
//			long long feature = 0;
//			unsigned char center = im[yc][xc];
//
//			for (int y = uu; y <= dd; y++) {
//				for (int x = ll; x <= rr; x++) {
//					feature |= ((long long)(im[y][x] > center) << idx);
//					idx++;
//				}
//			}
//
//			census[yc][xc] = feature;
//		}
//	}
//
//	return census;
//}
//
//VECBITMAP<float> ComputeCensusTensor(VECBITMAP<unsigned char>& imL, VECBITMAP<unsigned char>& imR, int sign)
//{
//	VECBITMAP<float> dsi(nrows, ncols, ndisps);
//	VECBITMAP<long long> censusL = ComputeCensusImage(imL);
//	VECBITMAP<long long> censusR = ComputeCensusImage(imR);
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			for (int d = 0; d < ndisps; d++) {
//				int xm = (x + sign*d + ncols) % ncols;
//				dsi.get(y, x)[d] = hamdist(censusL[y][x], censusR[y][xm]);
//			}
//		}
//	}
//	return dsi;
//}
//
//VECBITMAP<float> ComputeAdcensusTensor(cv::Mat& cvimL, cv::Mat& cvimR, int sign)
//{
//	const float ad_lambda = 30;
//	const float census_labmda = 10;
//
//	cv::Mat cvgrayL, cvgrayR;
//	cv::cvtColor(cvimL, cvgrayL, CV_BGR2GRAY);
//	cv::cvtColor(cvimR, cvgrayR, CV_BGR2GRAY);
//	VECBITMAP<unsigned char> grayL(nrows, ncols, 1, cvgrayL.data);
//	VECBITMAP<unsigned char> grayR(nrows, ncols, 1, cvgrayR.data);
//	VECBITMAP<float> dsi_census = ComputeCensusTensor(grayL, grayR, sign);
//
//	assert(cvimL.isContinuous());
//	assert(cvimR.isContinuous());
//	VECBITMAP<unsigned char> imL(nrows, ncols, 3, cvimL.data);
//	VECBITMAP<unsigned char> imR(nrows, ncols, 3, cvimR.data);
//
//	VECBITMAP<float> dsi(nrows, ncols, ndisps);
//	#pragma omp parallel for
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			for (int d = 0; d < ndisps; d++) {
//				int xm = (x + sign*d + ncols) % ncols;
//				unsigned char *pL = imL.get(y, x);
//				unsigned char *pR = imR.get(y, xm);
//				dsi.get(y, x)[d] = (fabs((float)pL[0] - pR[0])
//								+ fabs((float)pL[1] - pR[1])
//								+ fabs((float)pL[2] - pR[2])) / 3.f;
//			}
//		}
//	}
//
//	#pragma omp parallel for
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			for (int d = 0; d < ndisps; d++) {
//				dsi.get(y, x)[d] = 2 - exp(-dsi.get(y, x)[d] / ad_lambda) - exp(-dsi_census.get(y, x)[d] / census_labmda);
//			}
//		}
//	}
//
//	return dsi;
//}
//
//VECBITMAP<float> ComputeAdgradientTensor(VECBITMAP<float>& colgradL, VECBITMAP<float>& colgradR, int sign)
//{
//	VECBITMAP<float> dsi(nrows, ncols, ndisps);
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			for (int d = 0; d < ndisps; d++) {
//				int xm = (x + sign*d + ncols) % ncols;
//				float *pL = colgradL.get(y, x);
//				float *pR = colgradR.get(y, xm);
//				float cost_col  = fabs(pL[0] - pR[0]) + fabs(pL[1] - pR[1]) + fabs(pL[2] - pR[2]);
//				float cost_grad = fabs(pL[3] - pR[3]) + fabs(pL[4] - pR[4]);
//				cost_col = std::min((float)tau_col, cost_col);
//				cost_grad = std::min((float)tau_grad, cost_grad);
//				dsi.get(y, x)[d] = (1 - alpha) * cost_col + alpha * cost_grad;
//			}
//		}
//	}
//	return dsi;
//}
//
//cv::Mat SegmentToSuperPixels(cv::Mat& imL)
//{
//	cv::Mat argb(nrows, ncols, CV_8UC4);
//	assert(argb.isContinuous());
//	VECBITMAP<int> labelmap(nrows, ncols);
//
//	int from_to[] = { -1, 0, 0, 3, 1, 2, 2, 1 };
//	cv::mixChannels(&imL, 1, &argb, 1, from_to, 4);
//
//	int width(ncols), height(nrows), numlabels(0);;
//	unsigned int* pbuff = (unsigned int*)argb.data;
//	int* klabels = NULL;
//
//	int		k = 200;	// Desired number of superpixels.
//	double	m = 20;		// Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
//	SLIC segment;
//	segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, width, height, klabels, numlabels, k, m);
//	segment.DrawContoursAroundSegments(pbuff, klabels, width, height, 0xff0000);
//	delete[] klabels;
//
//	cv::Mat contour(nrows, ncols, CV_8UC3);
//	int from_to2[] = { 1, 2, 2, 1, 3, 0 };
//	cv::mixChannels(&argb, 1, &contour, 1, from_to2, 3);
//
//	return contour;
//}
//
//void prepareRect(int yc, int xc, cv::Rect& patchRect, cv::Rect& canvasRect)
//{
//	int yU = std::max(0, yc - patch_r);
//	int xL = std::max(0, xc - patch_r);
//	int yD = std::min(nrows - 1, yc + patch_r);
//	int xR = std::min(ncols - 1, xc + patch_r);
//	canvasRect = cv::Rect(xL, yU, xR - xL + 1, yD - yU + 1);
//	// From canvas coordinate to patch coordinate
//	patchRect = cv::Rect(xL - xc + patch_r, yU - yc + patch_r, xR - xL + 1, yD - yU + 1);
//}
//
//void on_mouse(int event, int x, int y, int flags, void *param)
//{
//
//	CvPoint left, right;
//	cv::Mat tmp;
//
//	if (event == CV_EVENT_MOUSEMOVE)
//	{
//		if (x >= ncols && x < 2 * ncols && y >= 0 && y < 2 * nrows)
//		{
//			tmp = g_stereo.clone();
//			left.x = x; left.y = y;
//			right.x = x + ncols; 	right.y = y;
//			cv::line(tmp, left, right, cv::Scalar(255, 0, 0, 1));
//
//			left.x = x; left.y = y;
//			right.x = x - ncols; 	right.y = y;
//			cv::line(tmp, left, right, cv::Scalar(255, 0, 0, 1));
//
//			cv::imshow("disparity", tmp);
//		}
//	}
//	if (event == CV_EVENT_RBUTTONDOWN)
//	{
//		if (x >= ncols && x < 3 * ncols)
//		{
//			tmp = g_stereo.clone();
//			x %= ncols; y %= nrows;
//			float gtDisp = g_GT.at<cv::Vec3b>(y, x)[0];
//			float myDisp = g_mydisp.at<cv::Vec3b>(y, x)[0];
//			printf("(%d, %d)->GT: %.1f, MINE: %.1f    ", y, x, gtDisp / scale, myDisp / scale);
//			printf("a = %.5f  b = %.5f  c = %.5f\n", g_coeffsL[y][x].a, g_coeffsL[y][x].b, g_coeffsL[y][x].c);
//
//			VECBITMAP<float> weights = ComputeLocalPatchWeights(y, x, g_imL);
//			float wsum = 0;  for (int i = 0; i < patch_w * patch_w; i++) wsum += weights.data[i];
//			float cost = ComputePlaneCost(y, x, g_coeffsL[y][x], g_colgradL, g_colgradR, weights, -1) / wsum;
//			printf("plane cost: %f\n", cost);
//
//			cv::Mat outImg1 = tmp(cv::Rect(0, 0, ncols, nrows)),
//				outImg2 = tmp(cv::Rect(ncols, 0, ncols, nrows)),
//				outImg3 = tmp(cv::Rect(0, nrows, ncols, nrows)),
//				outImg4 = tmp(cv::Rect(ncols, nrows, ncols, nrows)),
//				outImg5 = tmp(cv::Rect(ncols * 2, 0, ncols, nrows)),
//				outImg6 = tmp(cv::Rect(ncols * 2, nrows, ncols, nrows));
//			cv::Rect patchRect, canvasRect;
//			prepareRect(y, x, patchRect, canvasRect);
//
//			cv::Mat weightImg(patch_w, patch_w, CV_8UC3);
//			for (int y = 0; y < patch_w; y++) {
//				for (int x = 0; x < patch_w; x++) {
//					unsigned char w = weights[y][x] * 255.f;
//					weightImg.at<cv::Vec3b>(y, x) = cv::Vec3b(w, w, w);
//				}
//			}
//			(weightImg(patchRect)).copyTo(outImg5(canvasRect));
//
//			Plane coeff = g_coeffsL[y][x];
//			int yc = y, xc = x;
//			cv::Mat planeImg(patch_w, patch_w, CV_8UC3);
//			for (int y = 0; y < patch_w; y++) {
//				for (int x = 0; x < patch_w; x++) {
//					unsigned char d = scale * coeff.ToDisparity(y - patch_r + yc, x - patch_r + xc);
//					planeImg.at<cv::Vec3b>(y, x) = cv::Vec3b(d, d, d);
//				}
//			}
//			(planeImg(patchRect)).copyTo(outImg2(canvasRect));
//			cv::imshow("disparity", tmp);
//		}
//	}
//}
//
//void EvaluateDisparity(VECBITMAP<float>& h_disp, float thresh)
//{
//	float count[3] = { 0, 0, 0 };
//	float badPixelRate[3] = { 0, 0, 0 };
//	cv::Mat disp(nrows, ncols, CV_8UC3), badOnALL(nrows, ncols, CV_8UC3), badOnOC(nrows, ncols, CV_8UC3), gray;
//	cv::cvtColor(g_L, gray, CV_BGR2GRAY);
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//
//			count[0] += (g_OC.at<cv::Vec3b>(y, x)[0]   == 255);
//			count[1] += (g_ALL.at<cv::Vec3b>(y, x)[0]  == 255);
//			count[2] += (g_DISC.at<cv::Vec3b>(y, x)[0] == 255);
//
//			unsigned char g = gray.at<unsigned char>(y, x);
//			badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(g, g, g);
//			badOnOC.at<cv::Vec3b>(y, x)  = cv::Vec3b(g, g, g);
//
//			unsigned char d = (unsigned int)(scale * h_disp[y][x] + 0.5);
//			disp.at<cv::Vec3b>(y, x) = cv::Vec3b(d, d, d);
//
//			float diff = abs((float)disp.at<cv::Vec3b>(y, x)[0] - (float)g_GT.at<cv::Vec3b>(y, x)[0]);
//			if (g_OC.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh) {
//				badPixelRate[0]++;
//				badOnOC.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 255);
//			}
//			if (g_ALL.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh) {
//				badPixelRate[1]++;
//				badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(0, 0, 255);
//			}
//			if (g_DISC.at<cv::Vec3b>(y, x)[0] == 255 && diff > scale * thresh)	{
//				badPixelRate[2]++;
//			}
//	
//			
//
//			if (g_OC.at<cv::Vec3b>(y, x)[0] == 0) {		//draw occlusion region
//				badOnOC.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 0, 0);
//			}
//			if (g_ALL.at<cv::Vec3b>(y, x)[0] == 0) {	//draw occlusion region
//				badOnALL.at<cv::Vec3b>(y, x) = cv::Vec3b(255, 0, 0);
//			}
//		}
//	}
//
//	for (int i = 0; i < 3; i++) {
//		badPixelRate[i] /= count[i];
//	}
//	
//	printf("badPixelRate: %f%%, %f%%, %f%%\n",
//		badPixelRate[0] * 100.0f, badPixelRate[1] * 100.0f, badPixelRate[2] * 100.0f);
//	if (1)
//	{
//		g_segments = SegmentToSuperPixels(g_L);
//
//		cv::Mat compareImg(nrows * 2, ncols* 3, CV_8UC3);
//		cv::Mat outImg1 = compareImg(cv::Rect(0, 0, ncols, nrows)),
//			outImg2 = compareImg(cv::Rect(ncols, 0, ncols, nrows)),
//			outImg3 = compareImg(cv::Rect(0, nrows, ncols, nrows)),
//			outImg4 = compareImg(cv::Rect(ncols, nrows, ncols, nrows)),
//			outImg5 = compareImg(cv::Rect(ncols * 2, 0, ncols, nrows)),
//			outImg6 = compareImg(cv::Rect(ncols * 2, nrows, ncols, nrows));
//		g_GT.copyTo(outImg1);	disp.copyTo(outImg2);
//		g_segments.copyTo(outImg3);	badOnOC.copyTo(outImg4);
//		g_L.copyTo(outImg5);	badOnALL.copyTo(outImg6);
//
//		//cv::cvtColor(disp, disp, CV_RGB2GRAY);
//		cv::imwrite(folders[folder_id] + "disparity.png", disp);
//		cv::imwrite("d:\\disparity.png", disp);
//		cv::imwrite("d:\\badOnALL.png", badOnALL);
//		cv::imwrite("d:\\badOnOC.png", badOnOC);
//		
//
//		g_stereo = compareImg;
//		g_mydisp = disp;
//
//		cv::imshow("disparity", compareImg);
//		cv::setMouseCallback("disparity", on_mouse);	
//		cv::waitKey(0);
//
//	}
//}
//
//void TestSuperPixels()
//{
//	cv::Mat imL = cv::imread("./teddy/im2.png");
//	//cv::imshow("imL", imL);
//	//cv::waitKey(0);
//	nrows = imL.rows;
//	ncols = imL.cols;
//	cv::Mat argb(nrows, ncols, CV_8UC4);
//	cv::Mat contour(nrows, ncols, CV_8UC4);
//	assert(argb.isContinuous());
//
//
//	int from_to[] = { -1, 0, 0, 3, 1, 2, 2, 1 };
//	cv::mixChannels(&imL, 1, &argb, 1, from_to, 4);
//
//	int width(ncols), height(nrows);
//	unsigned int* pbuff = (unsigned int*)argb.data;
//	unsigned int* contour_buff = (unsigned int*)contour.data;
//
//	for (int iter = 0; iter < 5; iter++) {
//		int k = 200;//Desired number of superpixels.
//		double m = 20;//Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
//		cv::Mat labelmap(nrows, ncols, CV_32SC1);
//		int* klabels = (int*)labelmap.data;
//		int numlabels(0);
//
//		SLIC segment;
//		segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, width, height, klabels, numlabels, k, m);
//		memcpy(contour_buff, pbuff, nrows*ncols * 3 * sizeof(unsigned char));
//		segment.DrawContoursAroundSegments(contour_buff, klabels, width, height, 0xff0000);
//		if (klabels) delete[] klabels;
//
//		int from_to2[] = { 1, 2, 2, 1, 3, 0 };
//		cv::mixChannels(&contour, 1, &contour, 1, from_to2, 3);
//		imshow("contour", contour);
//		cv::waitKey(0);
//	}
//	if (pbuff) delete[] pbuff;
//	//printf("number of labels: %d\n", numlabels);
//}
//
//void PrepareColgradImages(cv::Mat& img, VECBITMAP<float>& colgrad)
//{
//	int sobel_scale = 1, sobel_delta = 0;
//	Mat gray, grad_x, grad_y;
//	
//	cv::cvtColor(img, gray, CV_BGR2GRAY);
//	cv::Sobel(gray, grad_x, CV_32F, 1, 0, 3, sobel_scale, sobel_delta, BORDER_DEFAULT);
//	cv::Sobel(gray, grad_y, CV_32F, 0, 1, 3, sobel_scale, sobel_delta, BORDER_DEFAULT);
//	//cv::convertScaleAbs(grad_x, grad_x);
//	//cv::convertScaleAbs(grad_y, grad_y);
//	grad_x = grad_x / 8.f;
//	grad_y = grad_y / 8.f;
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			colgrad.get(y, x)[0] = img.at<cv::Vec3b>(y, x)[0];
//			colgrad.get(y, x)[1] = img.at<cv::Vec3b>(y, x)[1];
//			colgrad.get(y, x)[2] = img.at<cv::Vec3b>(y, x)[2];
//			colgrad.get(y, x)[3] = grad_x.at<float>(y, x);
//			colgrad.get(y, x)[4] = grad_y.at<float>(y, x);
//		}
//	}
//}
//
//void RandomPermute(std::vector<cv::Point2d>& pointList, int k)
//{
//	int n = pointList.size();
//	for (int i = 0; i < k; i++) {
//		int j = rand() % n;
//		std::swap(pointList[i], pointList[j]);
//	}
//}
//
//inline double ComputePlaneCost(Plane& coeff, VECBITMAP<float>& dsi, std::vector<cv::Point2d>& pointList)
//{
//	double cost = 0;
//	int regionSize = pointList.size();
//	for (int i = 0; i < regionSize; i++) {
//		int y = pointList[i].y;
//		int x = pointList[i].x;
//		int d = 0.5 + coeff.a * x + coeff.b * y + coeff.c;
//		if (d < 0 || d > dmax) {
//			cost += BAD_PLANE_PENALTY;
//		}
//		else {
//			cost += dsi.get(y, x)[d];
//		}
//	}
//	return cost;
//}
//
//void RansacEstimate(std::vector<cv::Point2d>& pointList, VECBITMAP<float>& dsi, VECBITMAP<float>& disp, VECBITMAP<Plane>& coeffs)
//{
//	const int MIN_SAMPLE_SIZE = 5;
//	const int regionSize = pointList.size();
//
//	if (regionSize < MIN_SAMPLE_SIZE) {
//		for (int i = 0; i < regionSize; i++) {
//			int y = pointList[i].y;
//			int x = pointList[i].x;
//			coeffs[y][x].RandomAssign();
//		}
//		return;
//	}
//
//	cv::Mat A(MIN_SAMPLE_SIZE, 3, CV_64F);
//	cv::Mat b(MIN_SAMPLE_SIZE, 1, CV_64F);
//	cv::Mat c(3, 1, CV_64F);
//
//	Plane bestcoeff;
//	double bestcost = DBL_MAX;
//
//	for (int retry = 0; retry < 200; retry++) {
//
//		RandomPermute(pointList, MIN_SAMPLE_SIZE);
//		for (int i = 0; i < MIN_SAMPLE_SIZE; i++) {
//			int y = pointList[i].y;
//			int x = pointList[i].x;
//			A.at<double>(i, 0) = x;
//			A.at<double>(i, 1) = y;
//			A.at<double>(i, 2) = 1;
//			b.at<double>(i, 0) = disp[y][x];
//		}
//
//		//printf("solving\n");
//		cv::solve(A, b, c, DECOMP_SVD);
//		//printf("finished.\n");
//		Plane coeff;
//		coeff.a = c.at<double>(0, 0);
//		coeff.b = c.at<double>(1, 0);
//		coeff.c = c.at<double>(2, 0);
//
//		double cost = ComputePlaneCost(coeff, dsi, pointList);
//		if (cost < bestcost) {
//			bestcost = cost;
//			bestcoeff = coeff;
//		}
//	}
//
//	for (int i = 0; i < regionSize; i++) {
//		int y = pointList[i].y;
//		int x = pointList[i].x;
//		coeffs[y][x] = bestcoeff;
//	}
//}
//
//void RansacPlanefit(cv::Mat& imL, VECBITMAP<float>& dsiL, VECBITMAP<float>& dispL, VECBITMAP<Plane>& coeffsL)
//{
//	WinnerTakesAll(dsiL, dispL);
//
//	cv::Mat argb(nrows, ncols, CV_8UC4);
//	assert(argb.isContinuous());
//	VECBITMAP<int> labelmap(nrows, ncols);
//
//	int from_to[] = { -1, 0, 0, 3, 1, 2, 2, 1 };
//	cv::mixChannels(&imL, 1, &argb, 1, from_to, 4);
//	int width(ncols), height(nrows), numlabels(0);;
//	unsigned int* pbuff = (unsigned int*)argb.data;
//	int* klabels = NULL;
//
//	int		k = 200;	// Desired number of superpixels.
//	double	m = 20;		// Compactness factor. use a value ranging from 10 to 40 depending on your needs. Default is 10
//	SLIC segment;
//	segment.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(pbuff, width, height, klabels, numlabels, k, m);
//	segment.DrawContoursAroundSegments(pbuff, klabels, width, height, 0xff0000);
//
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			labelmap[y][x] = klabels[y*ncols + x];
//		}
//	}
//	delete[] klabels;
//
//
//	printf("numlabels: %d\n", numlabels);
//	std::vector<int> regionSize(numlabels, 0);
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			regionSize[labelmap[y][x]]++;
//		}
//	}
//
//
//	std::vector<std::vector<cv::Point2d>> regionList(numlabels);
//	for (int i = 0; i < numlabels; i++) {
//		regionList[i].reserve(regionSize[i]);
//	}
//
//
//	for (int y = 0; y < nrows; y++) {
//		for (int x = 0; x < ncols; x++) {
//			int label = labelmap[y][x];
//			regionList[label].push_back(cv::Point2d(x, y));
//		}
//	}
//	
//	#pragma omp parallel for
//	for (int id = 0; id < numlabels; id++) {
//		//printf("estimating region %d\n", id);
//		RansacEstimate(regionList[id], dsiL, dispL, coeffsL);
//	}
//	int toc = clock();
//	
//	PlaneMapToDisparityMap(coeffsL, dispL);
//	g_coeffsL = coeffsL;
//}
//
//int main()
//{
//#ifndef USE_OPENMP
//	omp_set_num_threads(1);
//#endif
//
//	Mat cvimL = imread(folders[folder_id] + "im2.png");
//	Mat cvimR = imread(folders[folder_id] + "im6.png");
//	nrows = cvimL.rows;
//	ncols = cvimL.cols;
//
//	assert(cvimL.isContinuous());
//	assert(cvimR.isContinuous());
//	VECBITMAP<unsigned char> imL(nrows, ncols, 3, cvimL.data);
//	VECBITMAP<unsigned char> imR(nrows, ncols, 3, cvimR.data);
//	VECBITMAP<float> dispL(nrows, ncols);
//	VECBITMAP<float> dispR(nrows, ncols);
//
//	VECBITMAP<float> colgradL(nrows, ncols, 5);
//	VECBITMAP<float> colgradR(nrows, ncols, 5);
//	PrepareColgradImages(cvimL, colgradL);
//	PrepareColgradImages(cvimR, colgradR);
//
//
//
//#ifdef USE_SUBPIXEL_MATCHING_COST
//
//	colgradL.SaveToBinaryFile("colgradL.bin");
//	colgradR.SaveToBinaryFile("colgradR.bin");
//
//	RunPatchMatchStereoSubpixel(imL, imR, colgradL, colgradR, dispL, dispR);
//#else
//
//
//	//Mat cvgrayL, cvgrayR;
//	//cvtColor(cvimL, cvgrayL, CV_BGR2GRAY);
//	//cvtColor(cvimR, cvgrayR, CV_BGR2GRAY);
//	//VECBITMAP<unsigned char> grayL(nrows, ncols, 1, cvgrayL.data);
//	//VECBITMAP<unsigned char> grayR(nrows, ncols, 1, cvgrayR.data);
//	//VECBITMAP<float> dsiL = ComputeCensusTensor(grayL, grayR, -1);
//	//VECBITMAP<float> dsiR = ComputeCensusTensor(grayR, grayL, +1);	
//
//	//VECBITMAP<float> dsiL = ComputeAdcensusTensor(cvimL, cvimR, -1);
//	//VECBITMAP<float> dsiR = ComputeAdcensusTensor(cvimR, cvimL, +1);
//
//	VECBITMAP<float> dsiL = ComputeAdgradientTensor(colgradL, colgradR, -1);
//	VECBITMAP<float> dsiR = ComputeAdgradientTensor(colgradR, colgradL, +1);
//
//	RunPatchMatchStereo(imL, imR, dsiL, dsiR, dispL, dispR);
//
//	//VECBITMAP<Plane> coeffsL(nrows, ncols);
//	//RansacPlanefit(cvimL, dsiL, dispL, coeffsL);
//	/*WinnerTakesAll(dsiL, dispL);*/
//#endif
//
//
//
//	g_L = cv::imread(folders[folder_id] + "im2.png");
//	g_GT = cv::imread(folders[folder_id] + "disp2.png");
//	g_OC = cv::imread(folders[folder_id] + "nonocc.png");
//	g_ALL = cv::imread(folders[folder_id] + "all.png");
//	g_DISC = cv::imread(folders[folder_id] + "disc.png");
//	g_colgradL = colgradL;
//	g_colgradR = colgradR;
//	g_imL = imL;
//	g_imR = imR;
//
//
//	//VECBITMAP<Plane> coeffsL(nrows, ncols);
//	//coeffsL.LoadFromBinaryFile(folders[folder_id] + "coeffsL.bin");
//	//PlaneMapToDisparityMap(coeffsL, dispL);
//	//g_coeffsL = coeffsL;
//
//	EvaluateDisparity(dispL, 0.5);
//
//	return 0;
//}
//
//
//
//
//
//
