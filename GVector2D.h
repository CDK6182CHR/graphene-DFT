/*
按OOP风格封装的gsl_vector
*/
#pragma once
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include<iostream>

class GVector2D
{
	gsl_vector* v;
public:
	GVector2D(double x = 0, double y = 0);
	GVector2D(const GVector2D& g);
	GVector2D(GVector2D&& g);
	inline gsl_vector* gvec() { return v; }
	inline const gsl_vector* gvec()const { return v; }
	inline double x()const { return gsl_vector_get(v, 0); }
	inline double y()const { return gsl_vector_get(v, 1); }
	inline double abs()const { return sqrt(x() * x() + y() * y()); }
	inline double square()const { return x() * x() + y() * y(); }
	inline double dis(const GVector2D& v)const {
		double dx = x() - v.x(), dy = y() - v.y();
		return sqrt(dx * dx + dy * dy);
	}
	double operator*(const GVector2D& v)const;//点乘
	GVector2D operator*(double d)const;
	GVector2D& operator=(const GVector2D& g);
	GVector2D& operator=(GVector2D&& g);
	GVector2D operator+(const GVector2D& g)const;
	GVector2D operator-(const GVector2D& g)const;
	GVector2D operator-()const;
	GVector2D& operator+=(const GVector2D& g);
	GVector2D& operator-=(const GVector2D& g);
	inline bool operator==(const GVector2D& g)const {
		return x() == g.x() && y() == g.y();
	}
	inline bool operator!=(const GVector2D& g)const {
		return x() != g.x() || y() != g.y();
	}
	inline bool notZero()const {
		return x() || y();
	}
	~GVector2D();
};


std::ostream& operator<<(std::ostream& out, const GVector2D& v);

