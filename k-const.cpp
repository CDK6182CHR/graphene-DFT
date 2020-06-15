#include"k-const.h"
#include"base.h"
#include<functional>
#include"GVector2D.h"
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_linalg.h>
#include<iostream>
using namespace std;


gsl_matrix_complex
* psi_1s1 = gsl_matrix_complex_alloc(RCount, RCount),
* psi_1s2 = gsl_matrix_complex_alloc(RCount, RCount),
* psi_2s1 = gsl_matrix_complex_alloc(RCount, RCount),
* psi_2s2 = gsl_matrix_complex_alloc(RCount, RCount);

gsl_matrix_complex* psi_all[NInnerOrbit] = {
	psi_1s1,psi_1s2,
	//psi_2s1,psi_2s2
};

GComplex
* kpsi_1s1 = new GComplex[NSet],
* kpsi_1s2 = new GComplex[NSet],
* kpsi_2s1 = new GComplex[NSet],
* kpsi_2s2 = new GComplex[NSet];

GComplex* kpsi_all[NInnerOrbit] = {
	kpsi_1s1,
	kpsi_1s2,
	//kpsi_2s1,
	//kpsi_2s2
};

gsl_matrix_complex
* Vopw_1s1 = gsl_matrix_complex_alloc(RCount,RCount),
* Vopw_1s2 = gsl_matrix_complex_alloc(RCount,RCount);
gsl_matrix_complex* Vopw_all[NInnerOrbit] = {
	Vopw_1s1,Vopw_1s2
};

gsl_matrix_complex* S = gsl_matrix_complex_alloc(NSet, NSet),
* Sinv = gsl_matrix_complex_alloc(NSet, NSet);

void cal_k_consts(const GVector2D& k)
{
	GComplex A(0, 0);  //相因子求和
	//for (int i = 0; i < N; i++) {
	//	for (int j = 0; j < N; j++) {
	//		A += gsl_complex_exp(gsl_complex_rect(0, (Rls[i] - Rls[j]) * k));
	//	}
	//}
	//A *= sqrt(Omega) / RCount / RCount / N;
	//cout << "A=" << A << endl;
	cal_psi_c(k, psi_1s1, phi_1s, r1, A);
	cal_psi_c(k, psi_1s2, phi_1s, r2, A);
	//cal_psi_c(k, psi_2s1, phi_2s, r1);
	//cal_psi_c(k, psi_2s2, phi_2s, r2);
	cal_k_psi_c(k, kpsi_1s1, psi_1s1,phi_1s1,A);
	cal_k_psi_c(k, kpsi_1s2, psi_1s2,phi_1s2,A);
	//cal_k_psi_c(k, kpsi_2s1, phi_2s1);
	//cal_k_psi_c(k, kpsi_2s2, phi_2s2);
	construct_S(k);
	cal_Vopw_matrix(k);
}

//对于每个k，对实空间元胞求和，计算psi_c
void cal_psi_c(const GVector2D& k,gsl_matrix_complex* psic, 
	std::function<double(double)> phi, const GVector2D& r0,const GComplex& A)
{
	gsl_matrix_complex_set_zero(psic);
	for (int a1 = 0; a1 < RCount; a1++) {
		for (int a2 = 0; a2 < RCount; a2++) {
			double r = dis(r0, a1, a2);
			GComplex p(0);
			for (int i = 0; i < N; i++) {
				p += gsl_complex_mul_real(gsl_complex_exp(gsl_complex_rect(0, k * Rls[i])),
					phi(r));
			}
			gsl_matrix_complex_set(psic, a1, a2, p / sqrt(N));
		}
	}
}

//对每个Kh计算  < k+Kh | psi_c >  数值积分
GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, 
	const gsl_matrix_complex* psi,const gsl_matrix* phic,const GComplex& A)
{
	GComplex p(0);
	for (int a1 = 0; a1 < RCount; a1++) {
		for (int a2 = 0; a2 < RCount; a2++) {
			GVector2D r = directPos(a1, a2);
			p += GComplex(gsl_complex_exp(gsl_complex_rect(0, -(Kh * r)-(k*r)))) *
				gsl_matrix_get(phic, a1, a2);
		}
		//todo: 这个积分表达式还不完全确定。
	}
	//static const double a = sqrt(N * Omega) / RCount / RCount;
	static const double a = sqrt(Omega) / RCount / RCount;
	return p * a;
}

//暂定用比较粗糙的矩形法数值积分公式计算
void cal_k_psi_c(const GVector2D& k, GComplex* kpsi, const gsl_matrix_complex* psic, 
	const gsl_matrix* phic, const GComplex& A)
{
	for (int i = 0; i < NSet; i++) {
		kpsi[i] = cal_k_Kh_psi_c(k, Khs[i], psic,phic,A);
	}
}

//构建S及其逆矩阵
void construct_S(const GVector2D& k)
{
	gsl_matrix_complex_set_zero(S);
	gsl_matrix_complex_set_zero(Sinv);
	static int signum;
	static gsl_permutation* P = gsl_permutation_alloc(NSet);
	for (int i = 0; i < NSet; i++) {
		for (int j = i; j < NSet; j++) {
			GComplex s ((int)(i == j),0);//矩阵元。首先有个Kroneck-delta
			for (int c = 0; c < NInnerOrbit; c++) {
				GComplex* kpsi_c = kpsi_all[c];
				s -= kpsi_c[i] * kpsi_c[j].conj();
			}
			gsl_matrix_complex_set(S, i, j, s);
			if (i != j) {
				gsl_matrix_complex_set(S, j, i, s.conj());  //Hermition性质
			}
		}
	}
	gsl_linalg_complex_LU_decomp(S, P, &signum);
	gsl_linalg_complex_LU_invert(S, P, Sinv);
}

void cal_Vopw_matrix(const GVector2D& k)
{
	for (int c = 0; c < NInnerOrbit; c++) {
		gsl_matrix_complex* vopw = Vopw_all[c];
		for (int i = 0; i < RCount; i++)
			for (int j = 0; j < RCount; j++)
				gsl_matrix_complex_set(vopw, i, j, Vopw(k, c, i, j));
	}
}

GComplex Vopw(const GVector2D& k,int c, int i, int j)
{
	double Ec = Ec_all[c];
	GComplex* kpsi_c = kpsi_all[c];
	const gsl_matrix_complex* psi_c = psi_all[c];
	GComplex res(0);
	for (int m = 0; m < NSet; m++) {
		const GVector2D& Kh = Khs[m];
		GVector2D&& r = directPos(i, j);
		res += kpsi_c[m] * gsl_complex_exp(gsl_complex_rect(0, (k * r) + (Kh * r))) *
			gsl_complex_conjugate(gsl_matrix_complex_get(psi_c, i, j));
	}
	static const double a = Omega*N;
	return res * Ec / a;
}




