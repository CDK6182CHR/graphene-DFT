#include "base.h"
#include<gsl/gsl_math.h>
#include<gsl/gsl_const.h>
#include<iostream>
#include<functional>
#include<stdint.h>
#include"GComplex.h"
using namespace std;

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
const double A0z = 1;  //z方向虚拟的晶格常数，归一化时候用
const double Omega = sqrt(3) * A0 * A0 / 2;

const double E1s = -me * pow(6 * e2k, 2) / (2 * hbar * hbar * 0.25);  //二维能量本征值
//const double E2s = -3.4 * GSL_CONST_MKSA_ELECTRON_CHARGE;
const double Ec_all[NInnerOrbit] = { 
	E1s,E1s,
	//E2s,E2s 
};
const int NValence = 8;//价电子数目
const double e = GSL_CONST_MKSA_ELECTRON_CHARGE;
const double e2k = e * e / (4 * M_PI * GSL_CONST_MKSA_VACUUM_PERMITTIVITY);  //e^2/(4 pi epsilon)

const int Z = 6;

inline int _cal_N_set();
//计算参数
//const double RCut=30e-10;//正空间晶格范围，即做FT的积分范围
const double KCut=8e10;//平面波截断半径，决定基组数目
const double prec = 1e-11; //收敛相对误差判据
const int MaxStep = 100;  //最大迭代步数
const int LHalfCount = 10;
const int LCount = 2 * LHalfCount + 1;
const int N = LCount * LCount;//晶胞数量
const int KCount=16;//1BZ高对称点路径每段折线的K点数目
const int RCount=40;//正空间元胞划分mesh的密度。将每一条基矢等分成多少段。
const int NSet = _cal_N_set();//基组数目


const GVector2D neigh[NNeigh] = {
	A1,-A1,A2,-A2,A1 + A2,A1 - A2,A2 - A1,-A1 - A2
};

const GVector2D KPath[KPOINTS] = {
	GVector2D(0,4 * M_PI / 3 / A0),  //K
	GVector2D(0,0),  //Gamma
	GVector2D(M_PI / sqrt(3) / A0,M_PI / A0),  //M
	GVector2D(0,4 * M_PI / 3 / A0),  //K
	//GVector2D(0,0),  //Gamma
};

GVector2D* Khs;

//inline int _cal_N() {
//	const int NX = RCut / A0;//最大的查找范围
//	int cnt = 0;
//	double A11 = A1.x(), A12 = A1.y(),
//		A21 = A2.x(), A22 = A2.y();
//	for (int i = -NX; i <= NX; i++)
//		for (int j = -NX; j <= NX; j++) {
//			double X = i * A11 + j * A21;
//			double Y = i * A12 + j * A22;
//			if (sqrt(X * X + Y * Y) <= RCut)
//				cnt++;
//		}
//	Rls = new GVector2D[cnt];
//	int t = 0;
//	for (int i = -NX; i <= NX; i++)
//		for (int j = -NX; j <= NX; j++) {
//			double X = i * A11 + j * A21;
//			double Y = i * A12 + j * A22;
//			if (sqrt(X * X + Y * Y) <= RCut) {
//				Rls[t++] = GVector2D(X, Y);
//			}
//		}
//	return cnt;
//}


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
			const GVector2D&& r = directPos(a, b);
			double s = phi(r.dis(r0));
			for (int i = 0; i < NNeigh; i++) {
				s += phi(r.dis(r0 + neigh[i]));
			}
			gsl_matrix_set(m, a, b, s);
		}
	}
	return m;
}

const gsl_matrix
* phi_1s1 = _init_phi_table(phi_1s, r1),
* phi_1s2 = _init_phi_table(phi_1s, r2);
//* phi_2s1 = _init_phi_table(phi_2s, r1),
//* phi_2s2 = _init_phi_table(phi_2s, r2);

double phi_1s(double r)
{
	//static const double aa = pow(ABohr, -1.5) / sqrt(M_PI);
	//return aa * exp(-r / ABohr);
	//按二维的波函数
	static const double aa = 4.0 / ABohr / sqrt(2 * M_PI);
	return aa * exp(-2 * r / ABohr);
}

double phi_2s(double r)
{
	return pow(ABohr, -1.5) * (1 - 0.5 * r / ABohr) * exp(-0.5 * r / ABohr) / sqrt(8 * M_PI);
}

//分数坐标->直角坐标。(a,b)/RCount为分数坐标
GVector2D directPos(int a, int b)
{
	return GVector2D((a+0.5)/RCount * A1.x() + (b+0.5)/RCount * A2.x(), 
		(a+0.5)/RCount * A1.y() + (b+0.5)/RCount * A2.y());
}

//double dis(const GVector2D& center, int a1, int a2)
//{
//	double X = (A1.x() * a1 + A2.x() * a2)/RCount;
//	double Y = (A1.y() * a1 + A2.y() * a2)/RCount;
//	double dx = center.x() - X, dy = center.y() - Y;
//	GVector2D dr(dx, dy);
//	
//	if (dr * A2 > A0 * A0 / 2)
//		dr -= A2;
//	else if (dr * A2 < -A0 * A0 / 2)
//		dr += A2;
//	if (dr * A1 > A0 * A0 / 2)
//		dr -= A1;
//	else if (dr * A1 < -A0 * A0 / 2)
//		dr += A1;
//
//	return dr.abs();
//}

//重新实现，按照最小像力约定，最近邻8个中找个最近的
//效率比较低！
double dis(const GVector2D& center, int a1, int a2)
{
	GVector2D r = directPos(a1, a2);
	double q = r.dis(center);
	for (int i = 0; i < NNeigh; i++) {
		double nq = r.dis(center + neigh[i]);
		if (nq < q)
			q = nq;
	}
	return q;
}

double dis(int a1, int b1, int a2, int b2)
{
	double da = double(a1 - a2) / RCount;  //a基矢方向投影
	double db = double(b1 - b2) / RCount;  //b基矢方向投影
	if (da > 0.5)da -= 1;
	else if (da < -0.5)da += 1;
	if (db > 0.5)db -= 1;
	else if (db < -0.5)db += 1;
	double dx = da * A1.x() + db * A2.x(), dy = da * A1.y() + db * A2.y();
	return sqrt(dx * dx + dy * dy);
}

double simpleDis(const GVector2D& center, int a1, int a2)
{
	double x = ((a1 + 0.5) * A1.x() + (a2 + 0.5) * A2.x()) / RCount;
	double y = ((0.5 + a1) * A1.y() + (0.5 + a2) * A2.y()) / RCount;
	double dx = x - center.x();
	double dy = y - center.y();
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

void output_reciprocal_matrix(const gsl_matrix_complex* m, const char* filename)
{
	FSTREAM(sf, filename, ios::out);
	for (int i = 0; i < NSet; i++) {
		for (int j = 0; j < NSet; j++)
			sf << GComplex(gsl_matrix_complex_get(m, i, j)) << '\t';
		sf << endl;
	}
	sf.close();
}

//计算离子势能
//如果效率许可，就直接用矢量运算
gsl_matrix* _init_Vext(const GVector2D& r0)
{
	gsl_matrix* m = gsl_matrix_alloc(RCount, RCount);
	for (int i = 0; i < RCount; i++) {
		for (int j = 0; j < RCount; j++) {
			const GVector2D&& r = directPos(i, j);
			const GVector2D&& t = r0 - r;
			double s = 0;
			for (int l1 = -LHalfCount; l1 <= LHalfCount; l1++) {
				for (int l2 = -LHalfCount; l2 <= LHalfCount; l2++) {
					const GVector2D&& Rl = A1 * l1 + A2 * l2;
					s += 1.0 / r.dis(r0 + Rl);
				}
			}
			gsl_matrix_set(m, i, j, -Z * e2k * s);
		}
	}
	return m;
}

const gsl_matrix* Vext_1 = _init_Vext(r1), * Vext_2 = _init_Vext(r2);


//按差向量为基本，打表所有点对的求和常数
gsl_matrix* _init_Vee_sum()
{
	gsl_matrix* m = gsl_matrix_alloc(2 * RCount - 1, 2 * RCount - 1);
	for (int a = -RCount + 1; a < RCount; a++) {
		for (int b = -RCount + 1; b < RCount; b++) {
			const GVector2D&& t = directPos(a, b);
			double s = 0;
			for (int l1 = -LHalfCount; l1 <= LHalfCount; l1++) {
				for (int l2 = -LHalfCount; l2 <= LHalfCount; l2++) {
					const GVector2D&& Rl = A1 * l1 + A2 * l2;
					//全0是对应相同位置且同一个元胞，干脆就不算，当成0处理
					if(a||b||l1||l2)
						s += 1.0 / (t + Rl).abs();
				}
			}
			gsl_matrix_set(m, a + RCount - 1, b + RCount - 1,
				e2k * Omega * s / RCount / RCount / 2);
		}
	}
	return m;
}

const gsl_matrix* sum_ee = _init_Vee_sum();
