/*对复值（实变）函数数值积分的支持*/
#pragma once
#include "GComplex.h"
#include "GVector2D.h"
#include<functional>

//复化simpson积分公式，对复值实变函数应用
GComplex simpson_complex(std::function<GComplex(double)> f, double a, double b, int N);

std::function<GComplex(double)> phi1s_factory(const GVector2D& q);