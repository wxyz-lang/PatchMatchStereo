#include <cassert>
#include <cstring>
#include <algorithm>

static const int n = 3;

struct NMPoint {
	float cost;
	float data[n];
	NMPoint operator+(NMPoint& m) { NMPoint ret(*this); for (int i = 0; i < n; i++) { ret.data[i] += m.data[i]; } return ret; }
	NMPoint operator-(NMPoint& m) { NMPoint ret(*this); for (int i = 0; i < n; i++) { ret.data[i] -= m.data[i]; } return ret; }
	NMPoint operator*(float c)    { NMPoint ret(*this); for (int i = 0; i < n; i++) { ret.data[i] *= c; } return ret; }
	NMPoint operator/(float c)    { NMPoint ret(*this); for (int i = 0; i < n; i++) { ret.data[i] /= c; } return ret; }
	NMPoint& operator+=(NMPoint& m) { for (int i = 0; i < n; i++) { this->data[i] += m.data[i]; } return *this; }
};

inline void ReplaceVertex(NMPoint *vertices, int idx, NMPoint& item )
{
	vertices[idx] = item;
	while (idx > 0 && vertices[idx - 1].cost > vertices[idx].cost) {
		std::swap(vertices[idx - 1], vertices[idx]);
		idx--;
	}

	for (int i = 0; i < n; i++) {
		assert(vertices[i].cost < vertices[i + 1].cost);
	}
}

inline void ReorderVertexList(NMPoint* vertices, int nitems)
{
	// Use simple insertion sort
	for (int i = 0; i < nitems; i++) {
		for (int j = i + 1; j < nitems; j++) {
			if (vertices[j].cost < vertices[i].cost) {
				std::swap(vertices[i], vertices[j]);
			}
		}
	}
	for (int i = 0; i < n; i++) {
		assert(vertices[i].cost < vertices[i + 1].cost);
	}
}

inline NMPoint ComputeGeometricCenter(NMPoint *vertices, int nitems)
{
	NMPoint sum(vertices[0]);
	for (int i = 1; i < nitems; i++) {
		sum += vertices[i];
	}
	return sum / nitems;
}

int NelderMeadOptimize(float *x, float(*feval)(float*, int), int maxiters = 0)
{
	int	retCode			= -1;
	const float tol		= 1.f;		// tol = 1 suffice.
	const float alpha	= 1;
	const float gamma	= 2;
	const float rho		= -0.5;
	const float sigma	= 0.5;

	NMPoint xo, xr, xe, xc, xBest, xWorst;
	NMPoint vertices[n + 1];

	for (int i = 0; i < n + 1; i++) {
		memcpy(vertices[i].data, x + i * n, n * sizeof(float));
		vertices[i].cost = feval(vertices[i].data, n);
	}
	ReorderVertexList(vertices, n + 1);

	float cost_before = vertices[0].cost;

	// FIXME: add convergence criteria
	// FIXME: add bound constraints

	for (int iter = 0; iter < maxiters; iter++) {

		// The list is ordered now
		xBest  = vertices[0];
		xWorst = vertices[n];
		if (xWorst.cost - xBest.cost < tol) {
			retCode = 0;
			break;
		}
		xo = ComputeGeometricCenter(vertices, n);

		// Reflection
		xr = xo + (xo - xWorst) * alpha;
		xr.cost = feval(xr.data, n);
		if (xBest.cost <= xr.cost && xr.cost < xWorst.cost) {
			ReplaceVertex(vertices, n, xr);
			continue;
		}

		// Expansion
		if (xr.cost < xBest.cost) {
			xe = xo + (xo - xWorst) * gamma;
			xe.cost = feval(xe.data, n);
			if (xe.cost < xr.cost) {
				ReplaceVertex(vertices, n, xe);
			}
			else {
				ReplaceVertex(vertices, n, xr);
			}
			continue;
		}

		// Contraction
		xc = xo + (xo - xWorst) * rho;
		xc.cost = feval(xc.data, n);
		if (xc.cost < xWorst.cost) {
			ReplaceVertex(vertices, n, xc);
			continue;
		}

		// Reduction
		for (int i = 1; i < n + 1; i++) {
			vertices[i] = xBest + (vertices[i] - xBest) * sigma;
			vertices[i].cost = feval(vertices[i].data, n);
		}
		ReorderVertexList(vertices, n + 1);
	}

	float cost_after = vertices[0].cost;
	assert(cost_after <= cost_before);
	if (cost_after > cost_before) {
		printf("FUCKING BUG, FUCKING BUG\n");
	}

	/*printf("Nelder-Mead returning.\n");*/
	for (int i = 0; i < n + 1; i++) {
		memcpy(x + i * n, vertices[i].data, n * sizeof(float));
	}

	return retCode;
}
