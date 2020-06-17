#include "integration.h"
#include "base.h"
#include <gsl/gsl_math.h>
#include <iostream>
using namespace std;

GComplex simpson_complex(std::function<GComplex(double)> f, double a, double b, int N)
{
	double h = (b - a) / N;
	GComplex s(0);
	for (int k = 0; k < N; k++) {
		s += f(a + (k + 0.5) * h) * 4 + f(a + k * h + h) * 2;
	}
	s += f(a) - f(b);
	s *= h / 6;
	return s;
}

std::function<GComplex(double)> phi1s_factory(const GVector2D& q)
{
	double qx = q.x(), qy = q.y();
	return [=](double x) {
		GComplex T(-2 / ABohr, -qx * cos(x) - qy * cos(x));
		return (T * T).reciprocal();
	};
}
