#pragma once

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>
#include <stack>

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


int main()
{
	FILE *fin = fopen("C:/Users/v-cz/Desktop/rank.txt", "r");
	FILE *fout = fopen("C:/Users/v-cz/Desktop/rank_latex.txt", "w");

	char rank[100];
	char buf1[100], buf2[100], buf3[100], buf4[100], buf5[100], buf6[100];

	while (fscanf(fin, "%s", rank) != EOF) {
		fprintf(fout, "& & $%s$ ", rank);
		for (int i = 0; i < 4; i++) {
			fscanf(fin, "%s %s %s %s %s %s", buf1, buf2, buf3, buf4, buf5, buf6);
			fprintf(fout, " & & $%s_{%s}$ & $%s_{%s}$ & $%s_{%s}$ ", buf1, buf2, buf3, buf4, buf5, buf6);
		}
		fprintf(fout, "\\\\\n");
	}

	fclose(fin);
	fclose(fout);
}