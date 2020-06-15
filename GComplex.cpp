#include"GComplex.h"

const GComplex GComplex::I(0, 1);
GComplex::GComplex(double r, double i) :c(gsl_complex_rect(r, i))
{
}

GComplex::GComplex(const gsl_complex& g) :
	c(gsl_complex_rect(GSL_REAL(g), GSL_IMAG(g)))
{

}

GComplex::operator gsl_complex() const
{
	return gsl_complex_rect(real(),imag());
}

GComplex GComplex::operator*(double x) const
{
	return GComplex(real() * x, imag() * x);
}

GComplex GComplex::operator*(const GComplex& g) const
{
	return GComplex(real() * g.real() - imag() * g.imag(),
		real() * g.imag() + imag()*g.real());
}

GComplex& GComplex::operator*=(double x)
{
	c.dat[0] *= x;
	c.dat[1] *= x;
	return *this;
}

GComplex& GComplex::operator+=(GComplex&& a)
{
	c.dat[0] += a.real();
	c.dat[1] += a.imag();
	return *this;
}

GComplex& GComplex::operator+=(gsl_complex&& a)
{
	c.dat[0] += GSL_REAL(a);
	c.dat[1] += GSL_IMAG(a);
	return *this;
}

GComplex& GComplex::operator-=(GComplex&& a)
{
	c.dat[0] -= a.real();
	c.dat[1] -= a.imag();
	return *this;
}

GComplex GComplex::operator+(const GComplex& a) const
{
	return GComplex(real() + a.real(), imag() + a.imag());
}

GComplex GComplex::operator-(const GComplex& a)const
{
	return GComplex(real() - a.real(), imag() - a.imag());
}

GComplex GComplex::operator-() const
{
	return GComplex(-real(),-imag());
}

GComplex GComplex::operator/(double x) const
{
	return GComplex(real() / x, imag() / x);
}

std::ostream& operator<<(std::ostream& out, const GComplex& g)
{
	out << "GComplex(" << g.real() << ", " << g.imag() << ')';
	return out;
}
