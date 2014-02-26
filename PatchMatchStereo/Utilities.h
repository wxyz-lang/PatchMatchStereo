#pragma once

#include <cmath>
#include <cstdlib>
#include <algorithm>

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
	void RandomAssign(int y, int x, int dmax)
	{
		const int RAND_HALF = RAND_MAX / 2;
		float z = dmax * ((double)rand() / RAND_MAX);
		float nx = ((double)rand() - RAND_HALF) / RAND_HALF;
		float ny = ((double)rand() - RAND_HALF) / RAND_HALF;
		float nz = ((double)rand() - RAND_HALF) / RAND_HALF;
		float norm = std::max(0.01f, sqrt(nx*nx + ny*ny + nz*nz));
		nx /= norm;
		ny /= norm;
		nz /= norm;
		*this = Plane(nx, ny, nz, y, x, z);
	}
	Plane RandomSearch(int y, int x, float radius_z0, float radius_n, float dmax)
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
#ifndef FRONTAL_PARALLEL_ONLY
		return Plane(nx__, ny__, nz__, y, x, z);
#else
		z = (int)(z + 0.5);
		return Plane(0.f, 0.f, 1.f, y, x, z);
#endif
	}
	float ToDisparity(int y, int x) { return a * x + b * y + c; }
	void GetAbc(float *v) { v[0] = a; v[1] = b; v[2] = c; }
	void SetAbc(float *abc)
	{
		a = abc[0];
		b = abc[1];
		c = abc[2];
	}
};

template<class T>
class VECBITMAP {
public:
	T *data;
	int w, h, n;
	bool is_shared;
	T *get(int y, int x) { return &data[(y*w + x)*n]; }		/* Get patch (y, x). */
	T *line_n1(int y) { return &data[y*w]; }				/* Get line y assuming n=1. */
	VECBITMAP() { w = h = n = 0; data = NULL; }
	VECBITMAP(const VECBITMAP& obj)
	{
		// This constructor is very necessary for returning an object in a function,
		// in the case that the Name Return Value Optimization (NRVO) is turned off.
		w = obj.w; h = obj.h; n = obj.n; is_shared = obj.is_shared;
		if (is_shared) { data = obj.data; }
		else { data = new T[w*h*n]; memcpy(data, obj.data, w*h*n*sizeof(T)); }
	}
	VECBITMAP(int h_, int w_, int n_ = 1, T* data_ = NULL)
	{
		w = w_; h = h_; n = n_;
		if (!data_) { data = new T[w*h*n]; is_shared = false; }
		else	    { data = data_;        is_shared = true; }
	}
	VECBITMAP& operator=(const VECBITMAP& m)
	{
		// printf("= operator invoked.\n");
		// FIXME: it's not suggested to overload assignment operator, should declare a copyTo() function instead.
		// However, if the assignment operator is not overloaded, do not invoke it (e.g. a = b), it is dangerous.
		if (data) { delete[] data; }
		w = m.w; h = m.h; n = m.n; is_shared = m.is_shared;
		if (m.is_shared) { data = m.data; }
		else { data = new T[w*h*n]; memcpy(data, m.data, w*h*n*sizeof(T)); }
		return *this;
	}
	~VECBITMAP() { if (!is_shared) delete[] data; }
	T *operator[](int y) { return &data[y*w]; }

	void SaveToBinaryFile(std::string filename)
	{
		FILE *fid = fopen(filename.c_str(), "wb");
		assert(fid != NULL);
		fwrite(data, sizeof(T), w*h*n, fid);
		fclose(fid);
	}
	void LoadFromBinaryFile(std::string filename)
	{
		FILE *fid = fopen(filename.c_str(), "rb");
		assert(fid != NULL);
		fread(data, sizeof(T), w*h*n, fid);
		fclose(fid);
	}
};


class Timer
{
public:
	static void tic()
	{
		time_stamps.push(clock());
	}
	static void tic(const char *msg)
	{
		printf("Processing %s ...\n", msg);
		tic();
	}
	static void toc()
	{
		clock_t tic = time_stamps.top();
		clock_t toc = clock();
		float time_elapsed = (toc - tic) / 1000.f;
		printf("%.2fs\n", time_elapsed);
		time_stamps.pop();
	}
private:
	static std::stack<clock_t> time_stamps;
};


void EvaluateDisparity(VECBITMAP<float>& h_disp, float thresh, VECBITMAP<Plane>& coeffsL = VECBITMAP<Plane>());
void RunLaplacianStereo(cv::Mat& imL, cv::Mat& imR, int ndisps);
VECBITMAP<float> ComputeAdGradientCostVolume(cv::Mat& imL, cv::Mat& imR, int ndisps, int sign, float granularity);
VECBITMAP<float> ComputeAdCensusCostVolume(cv::Mat& cvimL, cv::Mat& cvimR, int ndisps, int sign);
VECBITMAP<float> WinnerTakesAll(VECBITMAP<float>& dsi);

extern const std::string folders[60];
extern const int scale, ndisps, dmax, patch_w, patch_r, folder_id;
extern const float alpha, gamma, tau_col, tau_grad, granularity, BAD_PLANE_PENALTY;




