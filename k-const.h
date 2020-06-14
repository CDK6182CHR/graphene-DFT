/*包含所有对每个特定k点是常量的数据和函数*/
#pragma once
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<functional>
#include"GVector2D.h"
#include"GComplex.h"
#include "base.h"

extern gsl_matrix_complex * psi_1s1, * psi_1s2, * psi_2s1, * psi_2s2;//正空间元胞函数表
extern gsl_matrix_complex* psi_all[NInnerOrbit];
extern gsl_matrix_complex* S,*Sinv;//交叠矩阵 及其逆矩阵. 原则上外界只使用Sinv
extern GComplex* kpsi_1s1, * kpsi_1s2, *kpsi_2s1, *kpsi_2s2;//<k+Kh|psi_c> 下标为Kh
extern GComplex* kpsi_all[NInnerOrbit];

//完成对一个k点的所有常量计算
void cal_k_consts(const GVector2D& k);

void cal_psi_c(
	const GVector2D& k,
	gsl_matrix_complex* psic,
	std::function<double(double)> phi,
	const GVector2D& r0
	);

GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, const gsl_matrix* phic);

//计算k点对应的所有Kh的 <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,//存放结果
	const gsl_matrix* phic//phi_c 芯电子局域轨道
	);

void construct_S(const GVector2D& k);