#pragma once

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>


template<class T>
class VECBITMAP {
public:
	T *data;
	int w, h, n;
	bool is_shared;
	T *get(int y, int x) { return &data[(y*w + x)*n]; }		/* Get patch (y, x). */
	T *line_n1(int y) { return &data[y*w]; }				/* Get line y assuming n=1. */
	VECBITMAP() { data = NULL; }
	VECBITMAP(const VECBITMAP& obj)
	{
		printf("copy construction invoked.\n");
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
	//VECBITMAP& operator=(const VECBITMAP& m)
	//{
	//	printf("assignment opt invoked.\n");
	//	if (data) { delete[] data; }
	//	w = m.w; h = m.h; n = m.n; is_shared = m.is_shared;
	//	if (m.is_shared) { data = m.data; }
	//	else { data = new T[w*h*n]; memcpy(data, m.data, w*h*n*sizeof(T)); }
	//	return *this;
	//}
	VECBITMAP& assign(const VECBITMAP& m)
	{
		printf("assignment func invoked.\n");
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


VECBITMAP<float> testFunc(int nrows, int ncols)
{
	if (nrows <= 3){
		VECBITMAP<float> ret(nrows, ncols);
		printf("%x\n", ret.data);
		ret[0][0] = 3.1f;
		return ret;
	}
	else {
		VECBITMAP<float> ret(nrows + 1, ncols + 1);
		printf("%x\n", ret.data);
		ret[0][0] = 3.1f;
		return ret;
	}
}


int kkk = 23;


void test2(int kkk) {
	printf("%d\n", kkk);
	printf("%d\n", ::kkk);
}

int main()
{

	test2(1);

	////VECBITMAP<float> a(3, 4);
	////{
	////	VECBITMAP<float> b;
	////	b = a;
	////	printf("%x\n", a.data);
	////	printf("%x\n", b.data);
	////}
	////int val = 3;
	////int val(3);

	//VECBITMAP<float> m = testFunc(3, 4);
	//printf("%x\n", m.data);
	//

	////VECBITMAP<float> a = m;
	//////a.assign(m);
	////printf("%x\n", a.data);

	return 0;
}