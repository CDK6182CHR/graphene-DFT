#pragma once
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_math.h>
#include<iostream>

class GComplex
{
	gsl_complex c;
public:
	static const GComplex I;
	GComplex(double r = 0, double i = 0);
	GComplex(const gsl_complex& g);
	operator gsl_complex()const;
	GComplex operator*(double x)const;
	GComplex operator*(const GComplex& x)const;
	GComplex& operator*=(double x);
	GComplex& operator+=(GComplex&& a);
	GComplex& operator+=(gsl_complex&& a);
	GComplex& operator-=(const GComplex& a);
	GComplex operator+(const GComplex& a)const;
	GComplex operator-(const GComplex& a)const;
	GComplex operator-()const;
	GComplex operator/(double x)const;
	GComplex& operator=(const GComplex& g);
	inline double real()const { return GSL_REAL(c); }
	inline double imag()const { return GSL_IMAG(c); }
	inline double abs()const { return sqrt(abs_square()); }
	inline double abs_square()const { return real() * real() + imag() * imag(); }
	inline GComplex conj()const { return GComplex(real(), -imag()); }
	inline bool notZero()const { return real() || imag(); }
	//不可重载bool转换，否则造成多处【运算】有歧义。bool也是数值类型！
	GComplex reciprocal()const;//倒数
};

std::ostream& operator<<(std::ostream& out, const GComplex& g);