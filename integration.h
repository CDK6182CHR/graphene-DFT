/*�Ը�ֵ��ʵ�䣩������ֵ���ֵ�֧��*/
#pragma once
#include "GComplex.h"
#include "GVector2D.h"
#include<functional>

//����simpson���ֹ�ʽ���Ը�ֵʵ�亯��Ӧ��
GComplex simpson_complex(std::function<GComplex(double)> f, double a, double b, int N);

std::function<GComplex(double)> phi1s_factory(const GVector2D& q);