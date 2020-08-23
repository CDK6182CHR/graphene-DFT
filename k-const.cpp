#include"k-const.h"
#include"base.h"
#include<functional>
#include"GVector2D.h"
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_linalg.h>
#include "integration.h"
#include<iostream>
using namespace std;


gsl_matrix_complex
* psi_1s1 = nullptr,
* psi_1s2 = nullptr;
//* psi_2s1 = gsl_matrix_complex_alloc(RCount, RCount),
//* psi_2s2 = gsl_matrix_complex_alloc(RCount, RCount);

gsl_matrix_complex* psi_all[NInnerOrbit] = {
	psi_1s1,psi_1s2,
	//psi_2s1,psi_2s2
};

GComplex
* kpsi_1s1 = nullptr,
* kpsi_1s2 = nullptr,
* kpsi_2s1 = nullptr,
* kpsi_2s2 = nullptr;

GComplex* kpsi_all[NInnerOrbit] = {
	kpsi_1s1,
	kpsi_1s2,
	//kpsi_2s1,
	//kpsi_2s2
};

//gsl_matrix_complex
//* Vopw_1s1 = gsl_matrix_complex_alloc(RCount,RCount),
//* Vopw_1s2 = gsl_matrix_complex_alloc(RCount,RCount);
//gsl_matrix_complex* Vopw_all[NInnerOrbit] = {
//	Vopw_1s1,Vopw_1s2
//};

gsl_matrix_complex* S  ,
* Sinv;

void cal_k_consts(const GVector2D& k)
{
	cal_psi_c(k, psi_1s1,phi_1s1, phi_1s, r1);
	cal_psi_c(k, psi_1s2,phi_1s2, phi_1s, r2);
	cal_k_psi_c(k, kpsi_1s1, phi_1s, r1);
	cal_k_psi_c(k, kpsi_1s2, phi_1s, r2);
	construct_S(k);
	//cal_Vopw_matrix(k);
}

//对于每个k，对实空间元胞求和，计算psi_c
void cal_psi_c(const GVector2D& k,gsl_matrix_complex* psic, const gsl_matrix* phic,
	std::function<double(double)> phi, const GVector2D& r0)
{
	gsl_matrix_complex_set_zero(psic);
	for (int a1 = 0; a1< RCount; a1++) {
		for (int a2 = 0; a2 <RCount; a2++) {
			GComplex p(0);
			const GVector2D&& r = directPos(a1, a2);
			const GVector2D&& t = r0 - r;
			static const double aa = 4.0 / (ABohr * sqrt(2 * M_PI * N));
			for (int i = -LHalfCount; i <=LHalfCount; i++) {
				for (int j = -LHalfCount; j <=LHalfCount; j++) {
					const GVector2D&& Rl = A1 * i + A2 * j;
					//用解析形式的函数调用准确计算
					double d = (Rl + t).abs();
					p += GComplex(gsl_complex_exp(
						gsl_complex_rect(0,k * Rl)))*phi(d);
				}
			}
			gsl_matrix_complex_set(psic, a1, a2, p/LCount);
		}
	}
}

/*
旧版本实现的按实空间格点做0阶代数精度求和的算法，现已弃用
*/
# if 0
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
			//p += GComplex(gsl_complex_exp(gsl_complex_rect(0, -(k * r) - (Kh * r)))) *
			//	gsl_matrix_complex_get(psi, a1, a2);
		}
	}
	static const double a = sqrt(Omega) / RCount / RCount;
	return p * a;
}

//暂定用比较粗糙的矩形法数值积分公式计算
void cal_k_psi_c(const GVector2D& k, GComplex* kpsi, const gsl_matrix_complex* psic,
	const gsl_matrix* phic, const GComplex& A)
{
	for (int i = 0; i < NSet; i++) {
		kpsi[i] = cal_k_Kh_psi_c(k, Khs[i], psic, phic, A);
	}
}
#endif

GComplex cal_k_Kh_psi_c(const GVector2D& k, const GVector2D& Kh, 
	std::function<double(double)> phi, const GVector2D& r0)
{
	static const double aa = 4.0 / (ABohr * sqrt(2 * M_PI*Omega)) ;
	const GVector2D&& q = k + Kh;
	const auto&& f = phi1s_factory(k + Kh);
	static const int n = 5000;
	return simpson_complex(f, 0, 2 * M_PI, n) * aa * gsl_complex_exp(
		gsl_complex_rect(0, -q * r0)
		);
}




void cal_k_psi_c(const GVector2D& k, GComplex* kpsi, 
	const std::function<double(double)>& phi, const GVector2D& r0)
{
	for (int i = 0; i < NSet; i++) {
		kpsi[i] = cal_k_Kh_psi_c(k, Khs[i], phi, r0);
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
	//output_reciprocal_matrix(S, "S.txt");
	gsl_linalg_complex_LU_decomp(S, P, &signum);
	gsl_linalg_complex_LU_invert(S, P, Sinv);
}


/*
实空间坐标下等效的赝势值的计算。目前数量级上似存在问题
仅在需要在实空间算出等效赝势的时候才启用。目前选择直接计算H矩阵元中的赝势项。
*/
# if 0
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
	//cout << res * Ec  << endl;  
	return -res * Ec /sqrt(a);  
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

#endif


