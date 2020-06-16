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
//gsl_vector* Ek_sorted = gsl_vector_alloc(NSet);//�������������ֵ
gsl_matrix_complex* aKhs = gsl_matrix_complex_alloc(NSet, NSet);//����������һ��һ��

double V0;  //ƽ������

void init_density()
{
	gsl_matrix_complex_set_identity(aKhs);
	//for (int i = 0; i < RCount; i++)
	//	for (int j = 0; j < RCount; j++)
	//		gsl_matrix_set(density, i, j, NValence / Omega * RCount * RCount);
}

bool solve_k(const GVector2D& k)
{
	//init_density();  //for test
	cout << "Computing for k point " << k << endl;
	cal_k_consts(k);
	double Ecur = 0;  //��ǰ״̬����
	double En = 0;
	cout << "------------------------------------------" << endl;
	cout << "Step\t\ttotE" << endl;
	int step = 0;
	do {
		Ecur = En;//��һ�ε����������
		construct_Vr();
		construct_H(k);
		solve_single_secular(k);
		En = cal_density(k,step);
		cout << (step++) <<"\t\t"<<setprecision(14)<< En << endl;
		if (step >= MaxStep) {
			cout << "[Warning] Terminate forcely for k point " << k << endl;
			//break;  //todo: ���ԣ��������з�Χ����Ľ⣬�����������
			return false;
		}
	} while (fabs((En - Ecur) / En) > prec);
	cout << "converged for k point " << k << endl;
	//if (fabs(En) > 5e-17)
	//	return false;  //todo: ����ǿ��ȥ���쳣��
	return true;
}

//Ҫ���ǰVr�Ѿ������
void construct_H(const GVector2D& k)
{
	gsl_matrix_complex_set_zero(H);
	//�ȹ���Խ�Ԫ
	V0 = V_FT(GVector2D(0, 0)).real();
	cout << "\t\t\t\tV_0=" << V0 << endl;
	for (int i = 0; i < NSet; i++) {
		double h = (k + Khs[i]).square() * hbar_2me+V0;
		//cout << "V0=" << V0 << endl;
		//cout << "T+V0=" << h << endl;
		for (int c = 0; c < NInnerOrbit; c++) {
			const GComplex* kpsi = kpsi_all[c];
			h -= kpsi[i].abs_square() * Ec_all[c];
		}
		gsl_matrix_complex_set(H, i, i, gsl_complex_rect(h,0));
	}
	//�ǶԽ�Ԫ������Hermition����
	for (int i = 0; i < NSet; i++) {
		for (int j = i + 1; j < NSet; j++) {
			const GVector2D& Kh = Khs[i], Kh1 = Khs[j];
			GComplex h = V_FT(Kh - Kh1);
			//cout << "V_Kh=" << h << endl;
			for (int c = 0; c < NInnerOrbit; c++) {
				const GComplex* kpsi = kpsi_all[c];
				const GComplex& vps=kpsi[i] * kpsi[j].conj() * Ec_all[c];
				h -= vps;
				//cout << "V_PS=" << vps << endl;
			}
			
			gsl_matrix_complex_set(H, i, j, h);
			gsl_matrix_complex_set(H, j, i, h.conj());
		}
	}
}

//��ǰ�ܶ�Ӧ���Ѿ����
//�ݶ���0�״���������ֵ���ַ���
void construct_Vr()
{
	gsl_matrix_complex_set_zero(Vr);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GComplex v = 0;  //�õ��������ֵ
			double vee = 0;  //���ӻ���������
			GComplex vext = 0;  //�������ӻ���������
			const double p = gsl_matrix_get(density, i, j);
			const GVector2D r = directPos(i, j);
			//veeҪ����2�����ӻ�������Ϊһ�Ե���������
			vee = Vee(i, j);//���ӻ�������Ӧ�û�Ҫ���ϵ�ǰλ�õ��ܶ�
			//vext += VExtLocal(dis(r1, i, j),p);
			//vext += VExtLocal(dis(r2, i, j),p);
			vext += Veff1(i, j, 1);
			vext += Veff2(i, j, 1);
			//vext *= Omega / RCount / RCount;  //��û������ֵ���֣���Ҫ��һ����
			v = vext + VLDA(p) + vee;
			//cout << "vext=" << vext << endl << "VLDA=" << VLDA(p) << endl << "vee=" << vee << endl;
			//v = vext + VLDA(p);  //�����ӽ���
			//v = vext  + e2k * vee;  //Hatree����
			gsl_matrix_complex_set(Vr, i, j, v);
		}
	}
}

//�������ܡ����������Ԫ������
double Vext(const GVector2D& r0,int a,int b)
{
	static const double c0 = 7.4658;
	GVector2D&& r = directPos(a, b);
	double res = 1.0 / r.dis(r0);
	//for (int i = 0; i < NNeigh; i++) {
	//	const GVector2D& Rl = neigh[i];
	//	res += 1.0 / r.dis(r0 + Rl);
	//}
	res += c0 * LCount / A0;
	return res;
}

//��r1����Ϊ���ĵ���Ч�������ܼ��㣬�������Ƶĵֿ۲��֡�
inline GComplex Veff1(int a, int b,double p)
{
	double d = Vext(r1, a, b);
	return GComplex(- e2k * Z *d,0);
	//return -GComplex(gsl_matrix_complex_get(Vopw_1s1, a, b));
	//return GComplex(-N * e2k * Z / r) -gsl_matrix_complex_get(Vopw_1s1, a, b);
}

GComplex Veff1_debug(int a, int b, double p)
{
	double r = dis(r1, a, b);
	return GComplex(-N*e2k * Z / r,0);
	//return -GComplex(gsl_matrix_complex_get(Vopw_1s1, a, b));
	//return GComplex(- N*e2k * Z / r) - gsl_matrix_complex_get(Vopw_1s1, a, b);
}

inline GComplex Veff2(int a, int b,double p)
{
	double d = Vext(r2, a, b);
	return GComplex(- e2k * Z *d);
	//return GComplex(-p * e2k * Z / r) - gsl_matrix_complex_get(Vopw_1s2, a, b);
}

//���ӻ��������ܣ����������. ��������Vee
double Vee(int a, int b)
{
	double res = 0;
	const GVector2D&& r = directPos(a, b);
	static const double c0 = 7.4658;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			if (i != a || j != b) {
				const GVector2D&& r1 = directPos(i, j);  //r'
				double&& p = gsl_matrix_get(density, i, j);//����ʵ�ռ����ں���
				double v = 1.0 / r.dis(r1);
				//for (int l = 0; l < NNeigh; l++) {
				//	const GVector2D& Rl = neigh[l];
				//	v += 1.0 / r.dis(r1 + Rl);
				//}
				v += 7.4657 * LCount / A0;
				res += v * p;
			}
		}
	}
	return res * Omega / RCount / RCount * e2k / 2;
}

//������ֻ������ָ���ܴ��ĵ����ܶȣ�Ĭ�������򲢣��������۵��ӡ�
void cal_density_single(const GVector2D& k, int n)
{
	//static gsl_matrix* single_density = gsl_matrix_alloc(RCount, RCount);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GComplex psi(0);  //��������������ȡֵ
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
			double sq = psi.abs_square();
			gsl_matrix_set(density,i,j,sq+gsl_matrix_get(density,i,j));
		}
	}
}


//��չ��ϵ��ȷ���ռ����ܶ�
//��һ�ε���֮ǰ��aKhs�����ʼ������
//��ʱ�����򷽷��������
double cal_density(const GVector2D& k,int step)
{
	gsl_matrix_set_zero(density);
	gsl_eigen_hermv_sort(Ek, aKhs, GSL_EIGEN_SORT_VAL_ASC);
	static gsl_matrix* temp = gsl_matrix_alloc(NSet, NSet);
	if (step > 6) {//��ǰ�������ƽ��,��������
		for (int i = 0; i < NSet; i++)
			for (int j = 0; j < NSet; j++)
				gsl_matrix_set(temp, i, j, gsl_matrix_get(density, i, j));
	}
	double Etot = 0;
	//��֤�۵�����Ŀ��ż��
	for (int i = 0; i < NValence / 2; i++) {
		double Ei = gsl_vector_get(Ek, i);
		Etot += 2 * Ei;
		cal_density_single(k, i);
	}
	static const double rate = 1.00;
	static gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
	//unsigned long int Seed = 23410981;
	//gsl_rng_set(r, Seed);
	if (step > 6) {
		//double rate = gsl_rng_uniform(r) * 0.7 + 0.3;
		//double rate = 1.0;
		//cout << "rate=" << rate << endl;
		for (int i = 0; i < NSet; i++)
			for (int j = 0; j < NSet; j++)
				gsl_matrix_set(density, i, j,
					(rate*gsl_matrix_get(density, i, j) + (1-rate)*gsl_matrix_get(temp, i, j)));
	}
	norm_density();
	//////////////////
	double tot=0;
	for (int i = 0; i < RCount; i++)
		for (int j = 0; j < RCount; j++)
			tot += gsl_matrix_get(density, i, j);
	tot *= Omega / RCount / RCount;
	//cout << "totDensity: " << tot << endl;
	//////////////////
	return Etot;
}

//double VLDA(double p)
//{
//	//return 0;//todo ��ʱ���Ե����������ܣ�����Hatree����
//	static const double AA = -0.9164, BB = -0.2846, CC = 1.0529, DD = 0.3334,
//		FF = -0.096, GG = 0.0622, HH = -0.00232, JJ = 0.004;//CA-LDA����
//	double rs = pow(3 / (4 * M_PI * p), 1.0 / 3)*1e10;//����Angstrom��λ
//	double res = AA / rs;
//	if (rs >= 1) {
//		res += BB / (1 + CC * sqrt(rs) + DD * rs);
//	}
//	else {
//		res += FF + GG * log(rs) * HH * log(rs) + JJ * rs * log(rs);
//	}
//	double dif = -AA / rs / rs;  //���Ƕ�rs�ĵ�������
//	if (rs >= 1) {
//		dif -= BB / pow(1 + CC * sqrt(rs) + DD * rs, 2.0) * (CC / (2 * sqrt(rs) + DD));
//	}
//	else {
//		dif += GG / rs + HH + JJ * log(rs) + JJ;
//	}
//	res -= dif * pow(3.0 / (4 * M_PI), 1.0 / 3) * pow(p, -4.0 / 3) / 3.0;
//	return res*e;//���յ��ӷ���Ϊ��λ����
//}

double VLDA(double p)
{
	//ע��Ϊ���������⣬���˲�����
	static const double AA = -  e * e / (M_PI) * pow(3 * M_PI * M_PI, 1.0 / 2);
	return AA * sqrt(p);
}

GComplex V_FT(const GVector2D& Kh)
{
	GComplex res(0);
	
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GVector2D r = directPos(i, j);
			double p = gsl_matrix_get(density, i, j);
			GComplex v = gsl_matrix_complex_get(Vr, i, j);
			//if (Kh.notZero())
			//	v -= V0;//��ȥƽ������
			res += GComplex(v) * gsl_complex_exp(
				gsl_complex_rect(0, -(Kh * r)));
		}
	}
	//if (!Kh.notZero()) {
	//	cout << "Calculating V0 " << res << endl;
	//}
	return res / RCount/RCount;
}

//һ�������ڷ��̹��̵ķ�װ��S,H���Ѿ�����á�
void solve_single_secular(const GVector2D& k)
{
	static gsl_matrix_complex* C = gsl_matrix_complex_alloc(NSet, NSet);//�洢�м���
	//gsl_matrix_complex_set_zero(C);
	gsl_blas_zhemm(CblasLeft, CblasUpper, gsl_complex_rect(1, 0), Sinv, H,
		gsl_complex_rect(0, 0), C);
	//C=Sinv*H
	static gsl_eigen_hermv_workspace* ws = gsl_eigen_hermv_alloc(NSet);
	gsl_eigen_hermv(C, Ek, aKhs, ws);//�������������ҹ�һ����Ek�е�˳��һһ��Ӧ
}

void norm_density()
{
	double tot = 0;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			tot += gsl_matrix_get(density, i, j);
		}
	}
	tot *= Omega / RCount / RCount;
	//cout << "norm: tot=" << tot << endl;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			double p = gsl_matrix_get(density, i, j) / tot;
			gsl_matrix_set(density, i, j, NValence * p);
		}
	}
}
