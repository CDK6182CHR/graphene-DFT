#include "base.h"
#include<gsl/gsl_math.h>
#include<gsl/gsl_const.h>
#include<iostream>
#include<functional>
#include<stdint.h>
#include"GComplex.h"
using namespace std;

//inline gsl_vector* _vec_init(int flag) {
//	gsl_vector* v = gsl_vector_alloc(2);
//	switch (flag)
//	{
//	case 1:
//		gsl_vector_set(v, 0, sqrt(3) * A0 / 2);
//		gsl_vector_set(v, 1, A0 / 2); break;
//	case 2:
//		gsl_vector_set(v, 0, sqrt(3) * A0 / 2);
//		gsl_vector_set(v, 1, -A0 / 2); break;
//	case 3:
//		gsl_vector_set(v, 0, 2 * M_PI / A0 / sqrt(3));
//		gsl_vector_set(v, 1, 2 * M_PI / A0); break;
//	case 4:
//		gsl_vector_set(v, 0, 2 * M_PI / A0 / sqrt(3));
//		gsl_vector_set(v, 1, -2 * M_PI / A0); break;
//	case 5:
//		gsl_vector_set(v, 0, 0);
//		gsl_vector_set(v, 1, 0); break;
//	case 6:
//		gsl_vector_set(v, 0, A0 / sqrt(3));
//		gsl_vector_set(v, 1, 0); break;
//	default:
//		break;
//	}
//	return v;
//}

const GVector2D A1(sqrt(3)* A0 / 2, A0 / 2);
const GVector2D A2(sqrt(3)* A0 / 2, -A0 / 2);
const GVector2D B1(2 * M_PI / sqrt(3) / A0, 2 * M_PI / A0);
const GVector2D B2(2 * M_PI / sqrt(3) / A0, -2 * M_PI / A0);
const GVector2D r1(0, 0);
const GVector2D r2(A0 / sqrt(3), 0);

const double ABohr = GSL_CONST_MKSA_BOHR_RADIUS;  //Bohr半径0.529Angstrom
const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;  //约化Planck常数
const double me=GSL_CONST_MKSA_MASS_ELECTRON;//电子质量
const double hbar_2me = hbar * hbar / (2 * me);

const double A0=2.46e-10;//正空间晶格常数
const double Omega = sqrt(3) * A0 * A0 / 2;

const double E1s = -13.6 * GSL_CONST_MKSA_ELECTRON_CHARGE;
const double E2s = -3.4 * GSL_CONST_MKSA_ELECTRON_CHARGE;
const double Ec_all[NInnerOrbit] = { 
	E1s,E1s,
	//E2s,E2s 
};
const int NValence = 8;//价电子数目
const double e = GSL_CONST_MKSA_ELECTRON_CHARGE;
const double e2k = e * e / (4 * M_PI * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);  //e^2/(4 pi epsilon)

const int Z = 4;

inline int _cal_N();
inline int _cal_N_set();
//计算参数
const double RCut=30e-10;//正空间晶格范围，即做FT的积分范围
const double KCut=6e10;//平面波截断半径，决定基组数目
const int N=_cal_N();//晶胞数量
const int KCount=10;//1BZ高对称点路径每段折线的K点数目
const int RCount=40;//正空间元胞划分mesh的密度。将每一条基矢等分成多少段。
const int NSet = _cal_N_set();//基组数目

const GVector2D KPath[KPOINTS] = {
	GVector2D(0,4 * M_PI / 3 / A0),  //K
	GVector2D(0,0),  //Gamma
	GVector2D(M_PI / sqrt(3) / A0,M_PI / A0),  //M
	GVector2D(0,4 * M_PI / 3 / A0)  //K
};

GVector2D* Rls, * Khs;

inline int _cal_N() {
	const int NX = RCut / A0;//最大的查找范围
	int cnt = 0;
	double A11 = A1.x(), A12 = A1.y(),
		A21 = A2.x(), A22 = A2.y();
	for (int i = -NX; i <= NX; i++)
		for (int j = -NX; j <= NX; j++) {
			double X = i * A11 + j * A21;
			double Y = i * A12 + j * A22;
			if (sqrt(X * X + Y * Y) <= RCut)
				cnt++;
		}
	Rls = new GVector2D[cnt];
	int t = 0;
	for (int i = -NX; i <= NX; i++)
		for (int j = -NX; j <= NX; j++) {
			double X = i * A11 + j * A21;
			double Y = i * A12 + j * A22;
			if (sqrt(X * X + Y * Y) <= RCut) {
				Rls[t++] = GVector2D(X, Y);
			}
		}
	return cnt;
}

inline int _cal_N_set() {
	const int NX = KCut *0.5*A0/M_PI;//最大的查找范围
	cout << "Nx=" << NX << endl;
	int cnt = 0;
	double B11 = B1.x(), B12 = B1.y(),
		B21 = B2.x(), B22 = B2.y();
	for (int i = -NX; i <= NX; i++)
		for (int j = -NX; j <= NX; j++) {
			double X = i * B11 + j * B21;
			double Y = i * B12 + j * B22;
			if (sqrt(X * X + Y * Y) <= KCut)
				cnt++;
		}
	Khs = new GVector2D[cnt];
	int t = 0;
	for (int i = -NX; i <= NX; i++)
		for (int j = -NX; j <= NX; j++) {
			double X = i * B11 + j * B21;
			double Y = i * B12 + j * B22;
			if (sqrt(X * X + Y * Y) <= KCut) {
				Khs[t++] = GVector2D(X, Y);
			}
		}
	return cnt;
}

const gsl_matrix* _init_phi_table(function<double(double)> phi,
	const GVector2D& r0)
{
	gsl_matrix* m = gsl_matrix_alloc(RCount, RCount);
	for (int a = 0; a < RCount; a++) {
		for (int b = 0; b < RCount; b++) {
			double r = dis(r0, a, b);
			gsl_matrix_set(m, a, b, phi(r));
		}
	}
	return m;
}

const gsl_matrix
* phi_1s1 = _init_phi_table(phi_1s, r1),
* phi_1s2 = _init_phi_table(phi_1s, r2),
* phi_2s1 = _init_phi_table(phi_2s, r1),
* phi_2s2 = _init_phi_table(phi_2s, r2);

double phi_1s(double r)
{
	return pow(3 / ABohr, 1.5) * exp(-r / ABohr) / sqrt(M_PI);
	//return r;
}

double phi_2s(double r)
{
	return pow(ABohr, -1.5) * (1 - 0.5 * r / ABohr) * exp(-0.5 * r / ABohr) / sqrt(8 * M_PI);
}

GVector2D directPos(int a, int b)
{
	return GVector2D(a * A1.x() + b * A2.x(), a * A1.y() + b * A2.y());
}

double dis(const GVector2D& center, int a1, int a2)
{
	double X = (A1.x() * a1 + A2.x() * a2)/RCount;
	double Y = (A1.y() * a1 + A2.y() * a2)/RCount;
	double dx = center.x() - X, dy = center.y() - Y;
	GVector2D dr(dx, dy);
	
	if (dr * A2 > A0 * A0 / 2)
		dr -= A2;
	else if (dr * A2 < -A0 * A0 / 2)
		dr += A2;
	if (dr * A1 > A0 * A0 / 2)
		dr -= A1;
	else if (dr * A1 < -A0 * A0 / 2)
		dr += A1;

	return dr.abs();
}

double dis(int a1, int b1, int a2, int b2)
{
	double dx = double(a1 - a2) / RCount;
	double dy = double(b1 - b2) / RCount;
	return sqrt(dx * dx + dy * dy);
}

void output_real_matrix(const gsl_matrix* m, const char* filename)
{
	OPEN(fp, filename, "wb");
	int32_t r = RCount;
	fwrite(&r, sizeof(int32_t), 1, fp);
	fwrite(&r, sizeof(int32_t), 1, fp);
	for(int i=0;i<RCount;i++)
		for (int j = 0; j < RCount; j++) {
			double v = gsl_matrix_get(m, i, j);
			fwrite(&v, sizeof(double), 1, fp);
		}
	fclose(fp);
}

void output_real_matrix(const gsl_matrix_complex* m, const char* filename)
{
	OPEN(fp, filename, "wb");
	int32_t r = RCount;
	fwrite(&r, sizeof(int32_t), 1, fp);
	fwrite(&r, sizeof(int32_t), 1, fp);
	for (int i = 0; i < RCount; i++)
		for (int j = 0; j < RCount; j++) {
			const GComplex& v = gsl_matrix_complex_get(m, i, j);
			double d = v.abs();
			fwrite(&d, sizeof(double), 1, fp);
		}
	fclose(fp);
}
