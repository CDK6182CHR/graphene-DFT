/*�������ж�ÿ���ض�k���ǳ��������ݺͺ���*/
#pragma once
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<functional>
#include"GVector2D.h"
#include"GComplex.h"
#include "base.h"
#include<functional>

extern gsl_matrix_complex * psi_1s1, * psi_1s2, * psi_2s1, * psi_2s2;//���ռ�Ԫ��������
extern gsl_matrix_complex* psi_all[NInnerOrbit];
extern gsl_matrix_complex* S,*Sinv;//�������� ���������. ԭ�������ֻʹ��Sinv
extern GComplex* kpsi_1s1, * kpsi_1s2, *kpsi_2s1, *kpsi_2s2;//<k+Kh|psi_c> �±�ΪKh
extern GComplex* kpsi_all[NInnerOrbit];

extern gsl_matrix_complex* Vopw_1s1, * Vopw_1s2;
extern gsl_matrix_complex* Vopw_all[NInnerOrbit];

extern int NSet;//�������
extern GVector2D* Khs;  //����Ҫ���ǵĵ��ռ���

//��ɶ�һ��k������г�������
void cal_k_consts(const GVector2D& k);

void cal_psi_c(
	const GVector2D& k,
	gsl_matrix_complex* psic,
	const gsl_matrix* phic,  //���
	std::function<double(double)> phi,
	const GVector2D& r0
	);

//�ɰ汾ֱ����͵Ĵֲ�ʵ�֣���������
#if 0
GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, 
	const gsl_matrix_complex* psic, const gsl_matrix* phic,const GComplex& A);

//����k���Ӧ������Kh�� <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,//��Ž��
	const gsl_matrix_complex* psic,//psi_c о���Ӳ�����
	const gsl_matrix* phic,  //phi_c ����о���Ӳ�����
	const GComplex& A
	);
#endif

//������ֵ�����ֶ���
GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh,
	std::function<double(double)> phi, const GVector2D& r0);

//����ֵ���ַ�����������Kh�� <k+Kh | psi_c>
void cal_k_psi_c(
	const GVector2D& k,
	GComplex* kpsi,
	const std::function<double(double)>& phi,
	const GVector2D& r0
	);

void construct_S(const GVector2D& k);

//ʵ�ռ�����ֵ�Ķ������ݲ�����
#if 0
void cal_Vopw_matrix(const GVector2D& k);

GComplex Vopw(const GVector2D& k,int c/*о���ӹ�����*/, int i, int j);
#endif

int cal_N_set(const GVector2D& k);