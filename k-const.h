/*包含所有对每个特定k点是常量的数据和函数*/
#pragma once
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<functional>
#include"GVector2D.h"
#include"GComplex.h"
#include "base.h"
#include<functional>

extern gsl_matrix_complex * psi_1s1, * psi_1s2, * psi_2s1, * psi_2s2;//正空间元胞函数表
extern gsl_matrix_complex* psi_all[NInnerOrbit];
extern gsl_matrix_complex* S,*Sinv;//交叠矩阵 及其逆矩阵. 原则上外界只使用Sinv
extern GComplex* kpsi_1s1, * kpsi_1s2, *kpsi_2s1, *kpsi_2s2;//<k+Kh|psi_c> 下标为Kh
extern GComplex* kpsi_all[NInnerOrbit];

extern gsl_matrix_complex* Vopw_1s1, * Vopw_1s2;
extern gsl_matrix_complex* Vopw_all[NInnerOrbit];

extern int NSet;//基组个数
extern GVector2D* Khs;  //所有要考虑的倒空间格点

//完成对一个k点的所有常量计算
void cal_k_consts(const GVector2D& k);

void cal_psi_c(
	const GVector2D& k,
	gsl_matrix_complex* psic,
	const gsl_matrix* phic,  //打表
	std::function<double(double)> phi,
	const GVector2D& r0
	);

//旧版本直接求和的粗糙实现，现已弃用
#if 0
GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, 
	const gsl_matrix_complex* psic, const gsl_matrix* phic,const GComplex& A);

//计算k点对应的所有Kh的 <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,//存放结果
	const gsl_matrix_complex* psic,//psi_c 芯电子波函数
	const gsl_matrix* phic,  //phi_c 局域芯电子波函数
	const GComplex& A
	);
#endif

//运用数值积分手段算
GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh,
	std::function<double(double)> phi, const GVector2D& r0);

//用数值积分方法计算所有Kh的 <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,
	const std::function<double(double)>& phi,
	const GVector2D& r0
	);

void construct_S(const GVector2D& k);

//实空间赝势值的度量，暂不启用
#if 0
void cal_Vopw_matrix(const GVector2D& k);

GComplex Vopw(const GVector2D& k,int c/*芯电子轨道编号*/, int i, int j);
#endif

int cal_N_set(const GVector2D& k);