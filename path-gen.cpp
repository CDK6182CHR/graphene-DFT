#include "path-gen.h"
#include <iostream>
#include "base.h"
#include "k-const.h"
using namespace std;


//2020年8月23日新增系统初始化函数，负责空间分配之类的问题，
//防止分配前NSet之类的东西没初始化导致各种各样的问题。
void system_init()
{
	S = gsl_matrix_complex_alloc(NSet, NSet);
	Sinv = gsl_matrix_complex_alloc(NSet, NSet);
	kpsi_all[0] = kpsi_1s1 = new GComplex[NSet];
	kpsi_all[1] = kpsi_1s2 = new GComplex[NSet];
	psi_1s1 = gsl_matrix_complex_alloc(RCount, RCount);
	psi_1s2 = gsl_matrix_complex_alloc(RCount, RCount);
	init_Vtable();
}

void process_path()
{
	//采用C++面向对象的IO库，以便使用operator<<重载
	//OPEN(path, "KPATH.txt", "w");  //写入路径转折点表  （x,K空间坐标）
	//OPEN(fp, "EIGENVAL.txt", "w");  //横坐标点，所有能量本征值
	FSTREAM(path, "KPATH.txt", ios::out);
	FSTREAM(fp, "EIGENVAL.txt", ios::out);
	double x = 0;  //x轴位置
	int kcnt = 0;
	for (int i = 0; i < KPOINTS-1; i++) {
		const GVector2D& k1 = KPath[i], k2 = KPath[i + 1];
		cout << "KPATH-GEN: " << "k1=" << k1 << ", k2=" << k2 << endl;
		path << x << ' ' << k1 << endl;
		GVector2D k = k1;
		GVector2D dk = (k2 - k1) * (1.0/ double(KCount));
		for (int j = 0; j < KCount; j++) {
			cout << kcnt++ << ' ';
			bool f = solve_k(k);
			if (f) {
				fp << x << ' ';
				cout << "E[i]=";
				for (int h = 0; h < NSet; h++) {
					fp << gsl_vector_get(Ek, h) << ' ';
					//cout<< gsl_vector_get(Ek, h) << ' ';
				}
				fp << endl;
				fp.flush();
				cout << endl;
			}
			k += dk;
			x += 1;
		}
		if (i == KPOINTS - 2) {//最后一次循环，补充最后一个点
			cout << kcnt << ' ';
			bool f=solve_k(k);
			if (f) {
				fp << x << ' ';
				cout << "E[i]=";
				for (int h = 0; h < NSet; h++) {
					fp << gsl_vector_get(Ek, h) << ' ';
					//cout<< gsl_vector_get(Ek, h) << ' ';
				}
				fp << endl;
				cout << endl;
			}
		}
	}
	path << x << ' ' << KPath[KPOINTS - 1] << endl;
	path.close();
	fp.close();
}
