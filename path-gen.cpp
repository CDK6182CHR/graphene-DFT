#include "path-gen.h"
#include <iostream>
using namespace std;


void process_path()
{
	//OPEN(path, "KPATH.txt", "w");  //д��·��ת�۵��  ��x,K�ռ����꣩
	//OPEN(fp, "EIGENVAL.txt", "w");  //������㣬������������ֵ
	FSTREAM(path, "KPATH.txt", ios::out);
	FSTREAM(fp, "EIGENVAL.txt", ios::out);
	double x = 0;  //x��λ��
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
		if (i == KPOINTS - 2) {//���һ��ѭ�����������һ����
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
