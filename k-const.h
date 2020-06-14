/*�������ж�ÿ���ض�k���ǳ��������ݺͺ���*/
#pragma once
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<functional>
#include"GVector2D.h"
#include"GComplex.h"
#include "base.h"

extern gsl_matrix_complex * psi_1s1, * psi_1s2, * psi_2s1, * psi_2s2;//���ռ�Ԫ��������
extern gsl_matrix_complex* psi_all[NInnerOrbit];
extern gsl_matrix_complex* S,*Sinv;//�������� ���������. ԭ�������ֻʹ��Sinv
extern GComplex* kpsi_1s1, * kpsi_1s2, *kpsi_2s1, *kpsi_2s2;//<k+Kh|psi_c> �±�ΪKh
extern GComplex* kpsi_all[NInnerOrbit];

//��ɶ�һ��k������г�������
void cal_k_consts(const GVector2D& k);

void cal_psi_c(
	const GVector2D& k,
	gsl_matrix_complex* psic,
	std::function<double(double)> phi,
	const GVector2D& r0
	);

GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, const gsl_matrix* phic);

//����k���Ӧ������Kh�� <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,//��Ž��
	const gsl_matrix* phic//phi_c о���Ӿ�����
	);

void construct_S(const GVector2D& k);