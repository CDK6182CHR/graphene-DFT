#include "path-gen.h"
#include <iostream>
#include "base.h"
#include "k-const.h"
using namespace std;


//2020��8��23������ϵͳ��ʼ������������ռ����֮������⣬
//��ֹ����ǰNSet֮��Ķ���û��ʼ�����¸��ָ��������⡣
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
	//����C++��������IO�⣬�Ա�ʹ��operator<<����
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
		if (i == KPOINTS - 2) {//���һ��ѭ�����������һ����
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
