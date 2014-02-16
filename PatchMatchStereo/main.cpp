
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <vector>

#define PRECOMPUTE_WEIGHTS
#define DO_ITERATION
#define DO_POSTPROCESS


// Gobal variables
int nrows, ncols, ndisps, dmax;
float BAD_PLANE_PENALTY;
const int	patch_w				= 35;
const int	patch_r				= 17;
const int	maxiters			= 3;
const float gamma				= 10;
const float gamma_proximity		= 25;


int g_improve_cnt = 0;

struct Plane {
	float a, b, c;
	float nx, ny, nz;
	Plane() {}
	Plane(float a_, float b_, float c_, float nx_, float ny_, float nz_)
		:a(a_), b(b_), c(c_), nx(nx_), ny(ny_), nz(nz_) {}
	Plane(float nx_, float ny_, float nz_, int y, int x, float z)
	{
		const float eps = 0.001;
		nx = nx_;  ny = ny_;  nz = nz_;
		if (std::abs(nz) < eps) {
			if (nz > 0)  nz = +eps;
			else		 nz = -eps;
		}
		a = -nx / nz;
		b = -ny / nz;
		c = (nx * x + ny * y + nz * z) / nz;
	}
	Plane ReparametrizeInOtherView(int y, int x, int sign, int& qy, int &qx)
	{
		float z = a * x + b * y + c;
		qx = x + sign * z;
		qy = y;
		return Plane(nx, ny, nz, qy, qx, z);
	}
	Plane RandomSearch(int y, int x, float radius_z0, float radius_n)
	{
		const int RAND_HALF = RAND_MAX / 2;

		float nx__ = nx + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);
		float ny__ = ny + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);
		float nz__ = nz + radius_n * (((double)rand() - RAND_HALF) / RAND_HALF);

		float z = a * x + b * y + c
			+ radius_z0 * (((double)rand() - RAND_HALF) / RAND_HALF);
		z = std::max(0.f, std::min(z, (float)dmax));
		float norm = std::max(0.01f, sqrt(nx__*nx__ + ny__*ny__ + nz__*nz__));
		nx__ /= norm;
		ny__ /= norm;
		nz__ /= norm;

		return Plane(nx__, ny__, nz__, y, x, z);
	}
	float ToDisparity(int y, int x) { return a * x + b * y + c; }
};

template<class T>
class VECBITMAP {
public:
	T *data;
	int w, h, n;
	bool useOthersMemory;
	T *get(int y, int x) { return &data[(y*w + x)*n]; }		/* Get patch (x, y). */
	T *line_n1(int y) { return &data[y*w]; }				/* Get line y assuming n=1. */
	VECBITMAP() { }
	VECBITMAP(int h_, int w_, int n_ = 1, T* data_ = NULL)
	{
		w = w_; h = h_; n = n_;
		if (!data_) { data = new T[w*h*n]; useOthersMemory = false; }
		else	    { data = data_;        useOthersMemory = true; }
	}
	~VECBITMAP() { if (!useOthersMemory) delete[] data; }
	T *operator[](int y) { return &data[y*w]; }
};




inline bool InBound(int y, int x) { return 0 <= y && y < nrows && 0 <= x && x < ncols; }

float *PrecomputeWeights(VECBITMAP<unsigned char>& im)
{
	const int nw = nrows * ncols * patch_w * patch_w;
	const int patchsize = patch_w * patch_w;
	float *weights = new float[nw];
	memset(weights, 0, nw * sizeof(float));

	float w_prox[patch_w * patch_w];
	for (int y = 0; y < patch_w; y++) {
		for (int x = 0; x < patch_w; x++) {
			float dist = sqrt((y - patch_r) * (y - patch_r) + (x - patch_r) * (x - patch_r));
			w_prox[y * patch_w + x] = exp(-dist / gamma_proximity);
		}
	}

	for (int yc = 0; yc < nrows; yc++) {
		for (int xc = 0; xc < ncols; xc++) {

			//float *ww = &weights[(yc*ncols + xc) * patchsize];
			float *w = &weights[(yc*ncols + xc) * patchsize];
			unsigned char *rgb1 = im.get(yc, xc);
			unsigned char *rgb2 = NULL;

			int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
			int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

			for (int y = yb; y <= ye; y++) {
				//float *w = &ww[patch_w * (y - yc + patch_r)];
				for (int x = xb; x <= xe; x++) {
					rgb2 = im.get(y, x);
					float cost = std::abs((float)rgb1[0] - (float)rgb2[0])
							+ std::abs((float)rgb1[1] - (float)rgb2[1])
							+ std::abs((float)rgb1[2] - (float)rgb2[2]);
					//w[x - xc + patch_r] = cost;
					w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] = cost;
					//w[(y - yc + patch_r) * patch_w + (x - xc + patch_r)] 
					//	= exp(-cost / gamma) * w_prox[(y - yc + patch_r) * patch_w + (x - xc + patch_r)];
				}
			}
		}
	}

	for (int i = 0; i < nw; i++) {
		weights[i] = exp(-weights[i] / gamma);
	}

	return weights;
}

float ComputePlaneCost(int yc, int xc, Plane& coeff_try, VECBITMAP<float>& dsi, VECBITMAP<float>& w)
{
	float cost = 0.0f;
	for (int y = yc - patch_r; y <= yc + patch_r; y++) {
		for (int x = xc - patch_r; x <= xc + patch_r; x++) {
			int d = 0.5 + (coeff_try.a * x + coeff_try.b * y + coeff_try.c);
			if (InBound(y, x)) {
				if (d < 0 || d > dmax) {	// must be a bad plane.
					cost += BAD_PLANE_PENALTY;
				}
				else {
					cost += w[y - yc + patch_r][x - xc + patch_r] * (dsi.get(y, x)[d]);
				}
			}
		}
	}
	return cost;
}

void RandomInit(VECBITMAP<Plane>& coeffs, VECBITMAP<float>& bestcosts, VECBITMAP<float>& dsi, float* weights)
{
	const int	RAND_HALF = RAND_MAX / 2;
	const float eps = 0.001;

	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			// FIXME: Consider using double to avoid numerical issues
			float z = dmax * ((double)rand() / RAND_MAX);
			float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
			float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
			float nz = ((double)rand() - RAND_HALF) / RAND_HALF;

			float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
			nx /= norm;
			ny /= norm;
			nz /= norm;

			coeffs[y][x] = Plane(nx, ny, nz, y, x, z);
		}
	}

	int tic, toc;
	tic = clock();  printf("computing plane costs ...\n");	
	for (int y = 0; y < nrows; y++) {
		for (int x = 0; x < ncols; x++) {
			bestcosts[y][x] = ComputePlaneCost(y, x, coeffs[y][x], dsi, VECBITMAP<float>(patch_w, patch_w, 1, &weights[y*ncols + x]));
		}
	}
	toc = clock();  printf("compute plane costs comsumes %.2fs\n\n", (toc - tic) / 1000.f);
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

void ProcessView(VECBITMAP<Plane>& coeffsL, VECBITMAP<float>& bestcostsL, VECBITMAP<float>& dsiL, float* weightsL,
				 VECBITMAP<Plane>& coeffsR, VECBITMAP<float>& bestcostsR, VECBITMAP<float>& dsiR, float* weightsR,
	int iter, int sign)
{
	int xstart, xend, xchange, ystart, yend, ychange;
	if (iter % 2 == 0) {
		xstart = 0;  xend = ncols;  xchange = +1;
		ystart = 0;  yend = nrows;  ychange = +1;
	}
	else {
		xstart = ncols - 1;  xend = -1;  xchange = -1;
		ystart = nrows - 1;  yend = -1;  ychange = -1;
	}

	for (int y = ystart; y != yend; y += ychange) {
		for (int x = xstart; x != xend; x += xchange) {

			//VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[y*ncols + x]);
			VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y*ncols + x) * patch_w * patch_w]);

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
				Plane coeff_try = coeffsL[y][x].RandomSearch(y, x, radius_z, radius_n);
				ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
				radius_z /= 2.0f;
				radius_n /= 2.0f;
			}

			// View Propagation
			Plane coeff_try = coeffsL[y][x].ReparametrizeInOtherView(y, x, sign, qy, qx);
			if (0 <= qx && qx < ncols) {
				//VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[qy*ncols + qx]);
				VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(qy*ncols + qx) * patch_w * patch_w]);
				ImproveGuess(qy, qx, coeffsR[qy][qx], bestcostsR[qy][qx], coeff_try, dsiR, wR);
			}
		}
	}
}

VECBITMAP<float> ComputePatchWeights(int yc, int xc, VECBITMAP<unsigned char> img, float *weights)
{
	memset(weights, 0, patch_w * patch_w * sizeof(float));
	VECBITMAP<float> w(patch_w, patch_w, 1, weights);
	unsigned char *rgb1 = img.get(yc, xc);
	unsigned char *rgb2 = NULL;

	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

	for (int y = yb; y <= ye; y++) {
		for (int x = xb; x <= xe; x++) {
			rgb2 = img.get(y, x);
			float diff = std::abs((float)rgb1[0] - (float)rgb2[0])
					   + std::abs((float)rgb1[1] - (float)rgb2[1])
					   + std::abs((float)rgb1[2] - (float)rgb2[2]);
			w[y - yc + patch_r][x - xc + patch_r] = exp(-diff / gamma);
		}
	}

	return VECBITMAP<float>(patch_w, patch_w, 1, weights);
}

void ProcessView2(VECBITMAP<Plane>& coeffsL, VECBITMAP<float>& bestcostsL, VECBITMAP<float>& dsiL, VECBITMAP<unsigned char>& imL,
				  VECBITMAP<Plane>& coeffsR, VECBITMAP<float>& bestcostsR, VECBITMAP<float>& dsiR, VECBITMAP<unsigned char>& imR,
	int iter, int sign)
{
	int xstart, xend, xchange, ystart, yend, ychange;
	if (iter % 2 == 0) {
		xstart = 0;  xend = ncols;  xchange = +1;
		ystart = 0;  yend = nrows;  ychange = +1;
	}
	else {
		xstart = ncols - 1;  xend = -1;  xchange = -1;
		ystart = nrows - 1;  yend = -1;  ychange = -1;
	}

	float weights[patch_w *patch_w];

	for (int y = ystart; y != yend; y += ychange) {
		for (int x = xstart; x != xend; x += xchange) {

			VECBITMAP<float> wL = ComputePatchWeights(y, x, imL, weights);

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
				Plane coeff_try = coeffsL[y][x].RandomSearch(y, x, radius_z, radius_n);
				ImproveGuess(y, x, coeffsL[y][x], bestcostsL[y][x], coeff_try, dsiL, wL);
				radius_z /= 2.0f;
				radius_n /= 2.0f;
			}

			// View Propagation
			Plane coeff_try = coeffsL[y][x].ReparametrizeInOtherView(y, x, sign, qy, qx);
			if (0 <= qx && qx < ncols) {
				VECBITMAP<float> wR = ComputePatchWeights(qy, qx, imR, weights);
				ImproveGuess(qy, qx, coeffsR[qy][qx], bestcostsR[qy][qx], coeff_try, dsiR, wR);
			}
		}
	}
}

void WinnerTakesAll(VECBITMAP<float>& dsi, VECBITMAP<float> &disp)
{
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

void WeightedMedianFilter(int yc, int xc, VECBITMAP<float>& disp, VECBITMAP<float>& weights, VECBITMAP<bool>& valid)
{
	std::vector<std::pair<float, float>> dw_pairs;

	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

	for (int y = yb; y <= ye; y++) {
		for (int x = xb; x <= xe; x++) {
			if (valid[y][x]) {
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

void WeightedMedianFilter(int yc, int xc, VECBITMAP<float>& disp, VECBITMAP<float>& weights)
{
	std::vector<std::pair<float, float>> dw_pairs;

	int yb = std::max(0, yc - patch_r), ye = std::min(nrows - 1, yc + patch_r);
	int xb = std::max(0, xc - patch_r), xe = std::min(ncols - 1, xc + patch_r);

	for (int y = yb; y <= ye; y++) {
		for (int x = xb; x <= xe; x++) {
			std::pair<float, float> dw(disp[y][x], weights[y - yc + patch_r][x - xc + patch_r]);
			dw_pairs.push_back(dw);
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

void PostProcess(float *weightsL, VECBITMAP<Plane>& coeffsL, VECBITMAP<float>& dispL,
				 float *weightsR, VECBITMAP<Plane>& coeffsR, VECBITMAP<float>& dispR)
{
	//VECBITMAP<float> dispL(nrows, ncols), dispR(nrows, ncols);
	PlaneMapToDisparityMap(coeffsL, dispL);
	PlaneMapToDisparityMap(coeffsR, dispR);

	VECBITMAP<bool> validL(nrows, ncols), validR(nrows, ncols);
	CrossCheck(dispL, dispR, validL, validR);

#ifdef DO_POSTPROCESS
	// Hole filling
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

	PlaneMapToDisparityMap(coeffsL, dispL);
	PlaneMapToDisparityMap(coeffsR, dispR);

	for (int iter = 0; iter < 1; iter++) {

		CrossCheck(dispL, dispR, validL, validR);
		// Weighted median filtering
		for (int y = 0; y < nrows; y++) {
			if (y % 10 == 0) {
				printf("median filtering row %d\n", y);
			}
			for (int x = 0; x < ncols; x++) {
				if (!validL[y][x]){
					VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y * ncols + x) * patch_w * patch_w]);
					WeightedMedianFilter(y, x, dispL, wL);
				}
				if (!validR[y][x]){
					VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(y * ncols + x) * patch_w * patch_w]);
					WeightedMedianFilter(y, x, dispR, wR);
				}
			}
		}
	}

	CrossCheck(dispL, dispR, validL, validR);
	// Weighted median filtering
	for (int y = 0; y < nrows; y++) {
		if (y % 10 == 0) {
			printf("median filtering row %d\n", y);
		}
		for (int x = 0; x < ncols; x++) {
			if (!validL[y][x]){
				VECBITMAP<float> wL(patch_w, patch_w, 1, &weightsL[(y * ncols + x) * patch_w * patch_w]);
				WeightedMedianFilter(y, x, dispL, wL, validL);
			}
			if (!validR[y][x]){
				VECBITMAP<float> wR(patch_w, patch_w, 1, &weightsR[(y * ncols + x) * patch_w * patch_w]);
				WeightedMedianFilter(y, x, dispR, wR, validR);
			}
		}
	}
#endif
}

void RunPatchMatchStereo(VECBITMAP<unsigned char>& imL, VECBITMAP<float>& dsiL, VECBITMAP<Plane>& coeffsL, VECBITMAP<float>& dispL,
						 VECBITMAP<unsigned char>& imR, VECBITMAP<float>& dsiR, VECBITMAP<Plane>& coeffsR, VECBITMAP<float>& dispR)
{
	VECBITMAP<float> bestcostsL(nrows, ncols), bestcostsR(nrows, ncols);
	int tic, toc;
	FILE *fid = NULL;

#if defined(PRECOMPUTE_WEIGHTS)
	tic = clock();  printf("precomputing weights ...\n");
	float *weightsL = PrecomputeWeights(imL);
	float *weightsR = PrecomputeWeights(imR);
	toc = clock();  printf("precomputing weights use %.2fs\n\n", (toc - tic) / 1000.f);

	fid = fopen("weightsL.bin", "wb");
	fwrite(weightsL, sizeof(float), nrows * ncols * patch_w * patch_w, fid);
	fclose(fid);
	fid = fopen("weightsR.bin", "wb");
	fwrite(weightsR, sizeof(float), nrows * ncols * patch_w * patch_w, fid);
	fclose(fid);
#else
	const int nw = nrows * ncols * patch_w * patch_w;
	float *weightsL = new float[nw];
	float *weightsR = new float[nw];
	fid = fopen("weightsL.bin", "rb");
	fread(weightsL, sizeof(float), nrows * ncols * patch_w * patch_w, fid);
	fclose(fid);
	fid = fopen("weightsR.bin", "rb");
	fread(weightsR, sizeof(float), nrows * ncols * patch_w * patch_w, fid);
	fclose(fid); 
#endif


#if  defined(DO_ITERATION)
	// Random initialization
	tic = clock();  printf("performing random init ...\n");
	RandomInit(coeffsL, bestcostsL, dsiL, weightsL);
	RandomInit(coeffsR, bestcostsR, dsiR, weightsR);
	toc = clock();  printf("random init use %.2fs\n\n", (toc - tic) / 1000.f);

	// iteration
	for (int iter = 0; iter < maxiters; iter++) {

#ifdef PRECOMPUTE_WEIGHTS
		tic = clock();  printf("iter %d, scaning left view ...\n", iter);
		ProcessView(coeffsL, bestcostsL, dsiL, weightsL,
					coeffsR, bestcostsR, dsiR, weightsR,
					iter, -1);
		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);

		tic = clock();  printf("iter %d, scaning right view ...\n", iter);
		ProcessView(coeffsR, bestcostsR, dsiR, weightsR,
					coeffsL, bestcostsL, dsiL, weightsL,
					iter, +1);
		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
#else
		tic = clock();  printf("iter %d, scaning left view ...\n", iter);
		ProcessView2(coeffsL, bestcostsL, dsiL, imL,
					 coeffsR, bestcostsR, dsiR, imR,
					 iter, -1);
		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);

		tic = clock();  printf("iter %d, scaning right view ...\n", iter);
		ProcessView2(coeffsR, bestcostsR, dsiR, imR,
					 coeffsL, bestcostsL, dsiL, imL,
					 iter, +1);
		toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f);
#endif
	}
#endif

#ifdef DO_ITERATION
	fid = fopen("coeffsL.bin", "wb");
	fwrite(coeffsL.data, sizeof(Plane), nrows * ncols, fid);
	fclose(fid);
	fid = fopen("coeffsR.bin", "wb");
	fwrite(coeffsR.data, sizeof(Plane), nrows * ncols, fid);
	fclose(fid);
#else
	fid = fopen("coeffsL.bin", "rb");
	fread(coeffsL.data, sizeof(Plane), nrows * ncols, fid);
	fclose(fid);
	fid = fopen("coeffsR.bin", "rb");
	fread(coeffsR.data, sizeof(Plane), nrows * ncols, fid);
	fclose(fid);
#endif


	tic = clock();  printf("POST PROCESSING...\n");	
	PostProcess(weightsL, coeffsL, dispL, weightsR, coeffsR, dispR);
	toc = clock();  printf("%.2fs\n\n", (toc - tic) / 1000.f); 

	if (!weightsL) delete[] weightsL;
	if (!weightsR) delete[] weightsR;

	printf("g_improve_cnt: %d\n", g_improve_cnt);
}



int main()
{
	nrows = 375;
	ncols = 450;
	ndisps = 60;
	dmax = ndisps - 1;

	

	float *pdispL = new float[nrows * ncols];
	float *pdispR = new float[nrows * ncols];
	float *pdsiL = new float[nrows * ncols * ndisps];
	float *pdsiR = new float[nrows * ncols * ndisps];
	unsigned char *pimL = new unsigned char[nrows * ncols * 3];
	unsigned char *pimR = new unsigned char[nrows * ncols * 3];

	FILE *fid = NULL;
	fid = fopen("imL.bin", "rb");
	fread(pimL, sizeof(unsigned char), nrows * ncols * 3, fid);
	fclose(fid);
	fid = fopen("imR.bin", "rb");
	fread(pimR, sizeof(unsigned char), nrows * ncols * 3, fid);
	fclose(fid);
	fid = fopen("dsiL.bin", "rb");
	fread(pdsiL, sizeof(float), nrows * ncols * ndisps, fid);
	fclose(fid);
	fid = fopen("dsiR.bin", "rb");
	fread(pdsiR, sizeof(float), nrows * ncols * ndisps, fid);
	fclose(fid);

	int maxidx = 0, tensorsize = nrows * ncols * ndisps;
	for (int i = 0; i < ndisps; i++) {
		if (pdsiL[i] > pdsiL[maxidx]) {
			maxidx = i;
		}
	}
	//BAD_PLANE_PENALTY = 3 * pdsiL[maxidx];
	BAD_PLANE_PENALTY = 10;

	VECBITMAP<unsigned char> imL(nrows, ncols, 3, pimL);
	VECBITMAP<unsigned char> imR(nrows, ncols, 3, pimR);

	VECBITMAP<float> dispL(nrows, ncols, 1, pdispL);
	VECBITMAP<float> dispR(nrows, ncols, 1, pdispR);

	VECBITMAP<float> dsiL(nrows, ncols, ndisps, pdsiL);
	VECBITMAP<float> dsiR(nrows, ncols, ndisps, pdsiR);

	VECBITMAP<Plane> coeffsL(nrows, ncols);
	VECBITMAP<Plane> coeffsR(nrows, ncols);

	RunPatchMatchStereo(imL, dsiL, coeffsL, dispL, imR, dsiR, coeffsR, dispR);
	
	// Extract disparities
	//PlaneMapToDisparityMap(coeffsL, dispL);
	//PlaneMapToDisparityMap(coeffsR, dispR);



	fid = fopen("dispL.bin", "wb");
	fwrite(pdispL, sizeof(float), nrows * ncols, fid);
	fclose(fid);
	fid = fopen("dispR.bin", "wb");
	fwrite(pdispR, sizeof(float), nrows * ncols, fid);
	fclose(fid);

	delete[] pdispL, pdispR, pdsiL, pdsiR, pimL, pimR;
	return 0;
}
