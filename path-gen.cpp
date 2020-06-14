#include "path-gen.h"
#include <iostream>
using namespace std;


void process_path()
{
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
			solve_k(k);
			fp << x << ' ';
			cout << "E[i]=";
			for (int h = 0; h < NSet; h++) {
				fp << gsl_vector_get(Ek, h) << ' ';
				//cout<< gsl_vector_get(Ek, h) << ' ';
			}
			fp << endl;
			cout << endl;
			k += dk;
			x += dk.abs();
		}
		if (i == KPOINTS - 2) {//最后一次循环，补充最后一个点
			cout << kcnt << ' ';
			solve_k(k);
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
	path << x << ' ' << KPath[KPOINTS - 1] << endl;
	path.close();
	fp.close();
}
