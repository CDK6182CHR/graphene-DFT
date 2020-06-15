#include "GVector2D.h"
#include<gsl/gsl_complex_math.h>

GVector2D::GVector2D(double x, double y):
	v(gsl_vector_alloc(2))
{
	gsl_vector_set(v, 0, x);
	gsl_vector_set(v, 1, y);
}

GVector2D::GVector2D(const GVector2D& g):
	v(gsl_vector_alloc(2))
{
	gsl_vector_set(v, 0, gsl_vector_get(g.v,0));
	gsl_vector_set(v, 1, gsl_vector_get(g.v,1));
}

GVector2D::GVector2D(GVector2D&& g):
	v(g.v)
{
	g.v = nullptr;
}

double GVector2D::operator*(const GVector2D& v) const
{
	return x() * v.x() + y() * v.y();
}

GVector2D GVector2D::operator*(double d) const
{
	return GVector2D(x()*d,y()*d);
}

GVector2D& GVector2D::operator=(const GVector2D& g)
{
	gsl_vector_set(v, 0, g.x());
	gsl_vector_set(v, 1, g.y());
	return *this;
}

GVector2D& GVector2D::operator=(GVector2D&& g)
{
	v = g.v;
	g.v = nullptr;
	return *this;
}

GVector2D GVector2D::operator+(const GVector2D& g)const
{
	return GVector2D(x() + g.x(), y() + g.y());
}

GVector2D GVector2D::operator-(const GVector2D& g) const
{
	return GVector2D(x() - g.x(), y() - g.y());
}

GVector2D GVector2D::operator-() const
{
	return GVector2D(-x(),-y());
}

GVector2D& GVector2D::operator+=(const GVector2D& g)
{
	gsl_vector_set(v, 0, x() + g.x());
	gsl_vector_set(v, 1, y() + g.y());
	return *this;
}

GVector2D& GVector2D::operator-=(const GVector2D& g)
{
	gsl_vector_set(v, 0, x() - g.x());
	gsl_vector_set(v, 1, y() - g.y());
	return *this;
}

GVector2D::~GVector2D()
{
	if(v)
		gsl_vector_free(v);
	v = nullptr;
}

std::ostream& operator<<(std::ostream& out, const GVector2D& v)
{
	out << "GVector2D(" << v.x() << ", " << v.y() << ')';
	return out;
}
