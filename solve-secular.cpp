#include"solve-secular.h"
#include"base.h"
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_complex.h>
#include "k-const.h"
#include<gsl/gsl_math.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_rng.h>
#include<iostream>
#include<iomanip>
using namespace std;


gsl_matrix* density = gsl_matrix_alloc(RCount, RCount);
gsl_matrix_complex* Vr = gsl_matrix_complex_alloc(RCount, RCount);

gsl_matrix_complex* H = gsl_matrix_complex_alloc(NSet, NSet);

gsl_vector* Ek = gsl_vector_alloc(NSet);
//gsl_vector* Ek_sorted = gsl_vector_alloc(NSet);//经过排序的特征值
gsl_matrix_complex* aKhs = gsl_matrix_complex_alloc(NSet, NSet);//特征向量，一列一个

void init_density()
{
	gsl_matrix_complex_set_identity(aKhs);
}

bool solve_k(const GVector2D& k)
{
	//init_density();  //for test
	cout << "Computing for k point " << k << endl;
	cal_k_consts(k);
	construct_S(k);  //这是k点常量
	double Ecur = 0;  //当前状态能量
	double En = cal_density(k,0);
	cout << "------------------------------------------" << endl;
	cout << "Step\t\ttotE" << endl;
	int step = 0;
	do {
		Ecur = En;//上一次的总能量结果
		construct_Vr();
		construct_H(k);
		solve_single_secular(k);
		En = cal_density(k,step);
		cout << (step++) <<"\t\t"<<setprecision(14)<< En << endl;
		if (step >= MaxStep) {
			cout << "[Warning] Terminate forcely for k point " << k << endl;
			//break;  //todo: 测试，保留所有范围合理的解，无论收敛与否
			return false;
		}
	} while (fabs((En - Ecur) / En) > prec);
	cout << "converged for k point " << k << endl;
	//if (fabs(En) > 5e-17)
	//	return false;  //todo: 尝试强制去掉异常点
	return true;
}

//要求此前Vr已经构造好
void construct_H(const GVector2D& k)
{
	gsl_matrix_complex_set_zero(H);
	//先构造对角元
	double V0 = V_FT(GVector2D(0, 0)).real();
	cout << "\t\t\t\tV_0=" << V0 << endl;
	for (int i = 0; i < NSet; i++) {
		double h = (k + Khs[i]).square() * hbar_2me+V0;
		for (int c = 0; c < NInnerOrbit; c++) {
			const GComplex* kpsi = kpsi_all[c];
			//h -= kpsi[i].abs_square() * Ec_all[c];
		}
		gsl_matrix_complex_set(H, i, i, gsl_complex_rect(h,0));
	}
	//非对角元，按照Hermition构造
	for (int i = 0; i < NSet; i++) {
		for (int j = i + 1; j < NSet; j++) {
			const GVector2D& Kh = Khs[i], Kh1 = Khs[j];
			GComplex h = V_FT(Kh - Kh1);
			for (int c = 0; c < NInnerOrbit; c++) {
				const GComplex* kpsi = kpsi_all[c];
				//h -= kpsi[i] * kpsi[j].conj() * Ec_all[c];
			}
			gsl_matrix_complex_set(H, i, j, h);
			gsl_matrix_complex_set(H, j, i, h.conj());
		}
	}
}

//此前密度应该已经算好
//暂定用0阶代数精度数值积分方法
void construct_Vr()
{
	gsl_matrix_complex_set_zero(Vr);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GComplex v = 0;  //该点的势能数值
			double vee = 0;  //电子互作用势能
			GComplex vext = 0;  //电子离子互作用势能
			const double p = gsl_matrix_get(density, i, j);
			for (int a = 0; a < RCount; a++) {
				for (int b = 0; b < RCount; b++) {
					const double p1 = gsl_matrix_get(density, a, b);
					if (i != a || j != b) {
						vee += p1 / (dis(i, j, a, b));
					}
				}
			}
			vee *= p*Omega / (RCount * RCount-1);//电子互作用能应该还要算上当前位置的密度
			GVector2D r = directPos(i, j);
			//vext += VExtLocal(dis(r1, i, j),p);
			//vext += VExtLocal(dis(r2, i, j),p);
			vext += Veff1(i, j, p);
			vext += Veff2(i, j, p);
			//vext *= Omega / RCount / RCount;  //并没有做数值积分！不要归一化。
			v = vext + VLDA(p) + e2k * vee;
			gsl_matrix_complex_set(Vr, i, j, v);
		}
	}
}

//离子势能
double VExtLocal(double r,double p)
{
	static const double rc = A0 / 10;
	static const double AInner = 1.0 / rc;
	/*if (r > rc)
		return 1.0 / r;
	else
		return AInner;*/
	return -p*e2k * Z / r;
}

//以r1离子为中心的有效离子势能计算，包含赝势的抵扣部分。
GComplex Veff1(int a, int b,double p)
{
	double r = dis(r1, a, b);
	//return GComplex(-p * e2k * Z / r,0);
	//return -GComplex(gsl_matrix_complex_get(Vopw_1s1, a, b));
	return GComplex(-p * N*e2k * Z / r) - gsl_matrix_complex_get(Vopw_1s1, a, b);
}

GComplex Veff1_debug(int a, int b, double p)
{
	double r = dis(r1, a, b);
	return GComplex(-p * N*e2k * Z / r,0);
	//return -GComplex(gsl_matrix_complex_get(Vopw_1s1, a, b));
	return GComplex(-p * N*e2k * Z / r) - gsl_matrix_complex_get(Vopw_1s1, a, b);
}

GComplex Veff2(int a, int b,double p)
{
	double r = dis(r2, a, b);
	return GComplex(-p * e2k * Z / r) - gsl_matrix_complex_get(Vopw_1s2, a, b);
}

//本函数只叠加上指定能带的电子密度，默认自旋简并，即两个价电子。
void cal_density_single(const GVector2D& k, int n)
{
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GComplex psi(0);  //波函数在这个点的取值
			GVector2D r = directPos(i, j);
			for (int m = 0; m < NSet; m++) {
				GComplex a = gsl_matrix_complex_get(aKhs, m, n);
				const GVector2D& Kh = Khs[m];
				psi += a * gsl_complex_exp(gsl_complex_rect(0, k * r + Kh * r))
					/ sqrt(N * Omega);
				for (int c = 0; c < NInnerOrbit; c++) {
					GComplex* kpsi_c = kpsi_all[c];
					gsl_matrix_complex* psi_c = psi_all[c];
					psi -= (a * gsl_matrix_complex_get(psi_c, i, j)) *
						gsl_complex_conjugate(kpsi_c[m]);
				}
			}
			gsl_matrix_set(density, i, j, 2*psi.abs()+gsl_matrix_get(density,i,j));
		}
	}
}


//由展开系数确定空间电荷密度
//第一次调用之前，aKhs必须初始化！！
//暂时用排序方法解决问题
double cal_density(const GVector2D& k,int step)
{
	gsl_matrix_set_zero(density);
	gsl_eigen_hermv_sort(Ek, aKhs, GSL_EIGEN_SORT_VAL_ASC);
	static gsl_matrix* temp = gsl_matrix_alloc(NSet, NSet);
	if (step > 6) {//与前面的做个平均,加速收敛
		for (int i = 0; i < NSet; i++)
			for (int j = 0; j < NSet; j++)
				gsl_matrix_set(temp, i, j, gsl_matrix_get(density, i, j));
	}
	double Etot = 0;
	//保证价电子数目是偶数
	for (int i = 0; i < NValence / 2; i++) {
		double Ei = gsl_vector_get(Ek, i);
		Etot += 2 * Ei;
		cal_density_single(k, i);
	}
	//static const double rate = 0.50;
	static gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
	//unsigned long int Seed = 23410981;
	//gsl_rng_set(r, Seed);
	if (step > 6) {
		double rate = gsl_rng_uniform(r) * 0.7 + 0.3;
		//double rate = 1.0;
		//cout << "rate=" << rate << endl;
		for (int i = 0; i < NSet; i++)
			for (int j = 0; j < NSet; j++)
				gsl_matrix_set(density, i, j,
					(rate*gsl_matrix_get(density, i, j) + (1-rate)*gsl_matrix_get(temp, i, j)));
	}
	return Etot;
}

double VLDA(double p)
{
	//return 0;//todo 暂时忽略掉交换关联能，按照Hatree近似
	static const double AA = -0.9164, BB = -0.2846, CC = 1.0529, DD = 0.3334,
		FF = -0.096, GG = 0.0622, HH = -0.00232, JJ = 0.004;//CA-LDA常数
	double rs = pow(3 / (4 * M_PI * p), 1.0 / 3)*1e10;//按照Angstrom单位
	double res = AA / rs;
	if (rs >= 1) {
		res += BB / (1 + CC * sqrt(rs) + DD * rs);
	}
	else {
		res += FF + GG * log(rs) * HH * log(rs) + JJ * rs * log(rs);
	}
	double dif = -AA / rs / rs;  //这是对rs的导数部分
	if (rs >= 1) {
		dif -= BB / pow(1 + CC * sqrt(rs) + DD * rs, 2.0) * (CC / (2 * sqrt(rs) + DD));
	}
	else {
		dif += GG / rs + HH + JJ * log(rs) + JJ;
	}
	res -= dif * pow(3.0 / (4 * M_PI), 1.0 / 3) * pow(p, -4.0 / 3) / 3.0;
	return res*e;//按照电子伏特为单位处理
}

GComplex V_FT(const GVector2D& Kh)
{
	GComplex res(0);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GVector2D r = directPos(i, j);
			res += GComplex(gsl_matrix_complex_get(Vr, i, j)) * gsl_complex_exp(
				gsl_complex_rect(0,-(Kh * r)));
		}
	}
	return res / RCount/RCount;
}

//一次求解久期方程过程的封装。S,H都已经计算好。
void solve_single_secular(const GVector2D& k)
{
	static gsl_matrix_complex* C = gsl_matrix_complex_alloc(NSet, NSet);//存储中间结果
	//gsl_matrix_complex_set_zero(C);
	gsl_blas_zhemm(CblasLeft, CblasUpper, gsl_complex_rect(1, 0), Sinv, H,
		gsl_complex_rect(0, 0), C);
	//C=Sinv*H
	static gsl_eigen_hermv_workspace* ws = gsl_eigen_hermv_alloc(NSet);
	gsl_eigen_hermv(C, Ek, aKhs, ws);//特征向量正交且归一，和Ek中的顺序一一对应
}
