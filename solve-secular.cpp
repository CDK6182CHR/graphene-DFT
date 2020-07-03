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

gsl_matrix_complex* H = nullptr;

gsl_vector* Ek = nullptr;

gsl_matrix_complex* aKhs = nullptr;//����������һ��һ��


void init_density()
{
	//CA-LDA�Ե����ܶ��е����У������ܶȳ�����������׵���LDA��ը�����ʳ�ʼ��Ҫ��ɾ���
	//gsl_matrix_set_identity(density);
	for (int i = 0; i < RCount; i++)
		for (int j = 0; j < RCount; j++)
			gsl_matrix_set(density, i, j, NValence / Omega * RCount * RCount);
}

bool solve_k(const GVector2D& k)
{
	//init_density();  //���Գ�ʼ״̬��Ӱ��
	cout << "-------------------------------------------------------------" << endl;
	cout << "Computing for k point " << k << endl;
	cal_k_consts(k);
	double Ecur = 0;  //��ǰ״̬����
	double En = 0;
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
			return false;
		}
	} while (fabs((En - Ecur) / En) > prec);
	cout << "converged for k point " << k << endl;
	return true;
}

//Ҫ���ǰVr�Ѿ������
void construct_H(const GVector2D& k)
{
	gsl_matrix_complex_set_zero(H);
	//�ȹ���Խ�Ԫ
	double V0 = V_FT(GVector2D(0, 0)).real();
	cout << "\t\t\t\tV_0=" << V0 << endl;
	for (int i = 0; i < NSet; i++) {
		double h = (k + Khs[i]).square() * hbar_2me + V0;
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
			//veeҪ����2�����ӻ�������Ϊһ�Ե��������С����Ѿ�����ͳ��������֡�
			vee = Vee(i, j);
			vext += Veff1(i, j, 1);
			vext += Veff2(i, j, 1);
			//vext *= Omega / RCount / RCount;  //��û������ֵ���֣���Ҫ��һ����
			v = vext + VLDA(p) + vee;
			//cout << "vext=" << vext << endl << "VLDA=" << VLDA(p) << endl << "vee=" << vee << endl;
			gsl_matrix_complex_set(Vr, i, j, v);
		}
	}
}


//��r1����Ϊ���ĵ���Ч�������ܼ��㣬�������Ƶĵֿ۲��֡�
inline GComplex Veff1(int a, int b,double p)
{
	return GComplex(gsl_matrix_get(Vext_1, a, b));
}


inline GComplex Veff2(int a, int b,double p)
{
	return GComplex(gsl_matrix_get(Vext_2, a, b));
}

//���ӻ��������ܣ����������. ��������Vee
double Vee(int a, int b)
{
	double res = 0;
	const GVector2D&& r = directPos(a, b);
	static const double c0 = 7.4658;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			//if (i != a || j != b) {
			//	const GVector2D&& r1 = directPos(i, j);  //r'
			//	double&& p = gsl_matrix_get(density, i, j);//����ʵ�ռ����ں���
			//	//double v = dis(i,j, a, b);
			//	double v = 1.0 / r.dis(r1);
			//	//for (int l = 0; l < NNeigh; l++) {
			//	//	const GVector2D& Rl = neigh[l];
			//	//	v += 1.0 / r.dis(r1 + Rl);
			//	//}
			//	v += 7.4657 * LCount / A0;
			//	res += v * p;
			//}
			double&& p = gsl_matrix_get(density, i, j);
			res += p * Vee_sum_arg(i - a, j - b);
		}
	}
	return res * Omega / RCount / RCount * e2k / 2;
}

//������ֻ������ָ���ܴ��ĵ����ܶȣ�Ĭ�������򲢣��������۵��ӡ�
void cal_density_single(const GVector2D& k, int n)
{
	static gsl_matrix* single_density = gsl_matrix_alloc(RCount, RCount);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			GComplex psi(0);  //��������������ȡֵ
			GVector2D r = directPos(i, j);
			for (int m = 0; m < NSet; m++) {
				GComplex a = gsl_matrix_complex_get(aKhs, m, n);
				const GVector2D& Kh = Khs[m];
				psi += a * gsl_complex_exp(gsl_complex_rect(0, k * r + Kh * r))
					/ sqrt(Omega)/LCount;//todo: �˲���N?
				for (int c = 0; c < NInnerOrbit; c++) {
					GComplex* kpsi_c = kpsi_all[c];
					gsl_matrix_complex* psi_c = psi_all[c];
					psi -= (a * gsl_matrix_complex_get(psi_c, i, j)) *
						gsl_complex_conjugate(kpsi_c[m]);
				}
			}
			double sq = psi.abs_square();
			//gsl_matrix_set(density,i,j,sq+gsl_matrix_get(density,i,j));
			gsl_matrix_set(single_density, i, j, sq);
		}
	}
	norm_density_single(single_density, 2);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			gsl_matrix_set(density, i, j, gsl_matrix_get(density, i, j) +
				gsl_matrix_get(single_density, i, j));
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
	//static gsl_matrix* temp = gsl_matrix_alloc(NSet, NSet);
	//if (false&&(step+5)%10==0) {//��ǰ�������ƽ��,��������
	//	for (int i = 0; i < RCount; i++)
	//		for (int j = 0; j < RCount; j++)
	//			gsl_matrix_set(temp, i, j, gsl_matrix_get(density, i, j));
	//}
	double Etot = 0;
	//��֤�۵�����Ŀ��ż��
	for (int i = 0; i < NValence / 2; i++) {
		double Ei = gsl_vector_get(Ek, i);
		Etot += 2 * Ei;
		cal_density_single(k, i);
	}
	//norm_density();
	static const double rate = 1.0;
	static gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
	//unsigned long int Seed = 23410981;
	//gsl_rng_set(r, Seed);
	//if (false&&(step + 1) % 15 == 0) {
	//	//double rate = gsl_rng_uniform(r);
	//	//double rate = 1.0;
	//	cout << "Add average on density with rate=" << rate << endl;
	//	for (int i = 0; i < NSet; i++)
	//		for (int j = 0; j < NSet; j++)
	//			gsl_matrix_set(density, i, j,
	//				(rate*gsl_matrix_get(density, i, j) + 
	//					(1-rate)*gsl_matrix_get(temp, i, j)));
	//	//norm_density();
	//}
	
	//////////////////
	//double tot=0;
	//for (int i = 0; i < RCount; i++)
	//	for (int j = 0; j < RCount; j++)
	//		tot += gsl_matrix_get(density, i, j);
	//tot *= Omega / RCount / RCount;
	//cout << "totDensity: " << tot << endl;
	//////////////////
	return Etot;
}

double VLDA(double p)
{
	//return 0;//todo ��ʱ���Ե����������ܣ�����Hatree����
	if (p == 0)
		return 0;
	static const double AA = -0.9164, BB = -0.2846, CC = 1.0529, DD = 0.3334,
		FF = -0.096, GG = 0.0622, HH = -0.00232, JJ = 0.004;//CA-LDA����
	//double rs = pow(3 / (4 * M_PI * p), 1.0 / 3)*1e10;//����Angstrom��λ
	static const double KK = pow(9 * M_PI / 4, 1.0 / 3);
	double rs = KK / (sqrt(2 * M_PI * p))*1e10;
	double res = AA / rs;
	if (rs >= 1) {
		res += BB / (1 + CC * sqrt(rs) + DD * rs);
	}
	else {
		res += FF + GG * log(rs) * HH * log(rs) + JJ * rs * log(rs);
	}
	double dif = -AA / rs / rs;  //���Ƕ�rs�ĵ�������
	if (rs >= 1) {
		dif -= BB / pow(1 + CC * sqrt(rs) + DD * rs, 2.0) * (CC / (2 * sqrt(rs) + DD));
	}
	else {
		dif += GG / rs + HH + JJ * log(rs) + JJ;
	}
	res -= dif * pow(3.0 / (4 * M_PI), 1.0 / 3) * pow(p, -4.0 / 3) / 3.0;
	return res*e;//���յ��ӷ���Ϊ��λ����
}

//��һ����ʽ��LDA���������ƣ�����VASP-tutor
#if 0
double VLDA(double p)
{
	//ע��Ϊ���������⣬���˲�����
	//return 0;
	static const double AA = -  e * e / (M_PI) * pow(3 * M_PI * M_PI, 1.0 / 2);
	return AA * sqrt(p);
}
#endif

GComplex V_FT(const GVector2D& Kh)
{
	GComplex res(0);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			const GVector2D&& r = directPos(i, j);
			GComplex v = gsl_matrix_complex_get(Vr, i, j);
			res += GComplex(v) * gsl_complex_exp(
				gsl_complex_rect(0, -(Kh * r)));
		}
	}
	return res / RCount / RCount;
}
gsl_matrix_complex* C = nullptr;//�洢�м���

//һ�������ڷ��̹��̵ķ�װ��S,H���Ѿ�����á�
void solve_single_secular(const GVector2D& k)
{
	//gsl_matrix_complex_set_zero(C);
	gsl_blas_zhemm(CblasLeft, CblasUpper, gsl_complex_rect(1, 0), Sinv, H,
		gsl_complex_rect(0, 0), C);
	//C=Sinv*H
	gsl_eigen_hermv_workspace* ws = gsl_eigen_hermv_alloc(NSet);
	gsl_eigen_hermv(C, Ek, aKhs, ws);//�������������ҹ�һ����Ek�е�˳��һһ��Ӧ
	gsl_eigen_hermv_free(ws);
}

/*�������ܵ��ܶ��ٹ�һ���ǲ����ʵģ��������¸���������׿��ܲ����ȡ����á�*/
# if 0
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
#endif

void norm_density_single(gsl_matrix* m,int electrons)
{
	double tot = 0;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			tot += gsl_matrix_get(m, i, j);
		}
	}
	tot *= Omega / RCount / RCount;
	//cout << "norm: tot=" << tot << endl;
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			double p = gsl_matrix_get(m, i, j) / tot;
			gsl_matrix_set(m, i, j, electrons * p);
		}
	}
}
