#pragma once
/*���������ȵ�����*/
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex.h>
#include "GVector2D.h"
#include<iostream>
#include<fstream>

//��������
#define NInnerOrbit 2  //о���ӹ������
extern const GVector2D A1, A2, B1, B2;//���ռ�͵��ռ��ʸ
extern const GVector2D r1, r2;//����ԭ�ӵ�λ������
extern const double ABohr;  //Bohr�뾶0.529Angstrom
extern const double hbar;  //Լ��Planck����
extern const double me;//��������
extern const double hbar_2me;
extern const double Omega;//���ռ�Ԫ������������

extern const double A0;//���ռ侧����

extern const double E1s, E2s;  //ԭ�ӹ����������
extern const double Ec_all[NInnerOrbit];  //����ԭ�ӹ������������ֵ����������ʽ
extern const int NValence; //�۵���
extern const int Z;  //�������Ӻ˵������ֻ��۵���
extern const double e,e2k;  //Ԫ���


//�������
extern const double RCut;//���ռ侧��Χ������FT�Ļ��ַ�Χ
extern const double KCut;//ƽ�沨�ضϰ뾶
extern const int N;//��������
extern const int KCount;//1BZ�߶ԳƵ�·��ÿ�����ߵ�K����Ŀ
extern const int RCount;//���ռ仮��mesh���ܶȡ���ÿһ����ʸ�ȷֳɶ��ٶΡ�
extern const int NSet;//�������

extern GVector2D* Rls;  //����Ҫ���ǵ����ռ���
extern GVector2D* Khs;  //����Ҫ���ǵĵ��ռ���

extern const gsl_matrix* phi_1s1, * phi_1s2, * phi_2s1, * phi_2s2; //����о���Ӳ��������

#define KPOINTS 4  //���ռ�·�������



extern const GVector2D KPath[KPOINTS];

double phi_1s(double r);  //1s���������������ͬ��
double phi_2s(double r);  //2s���������������ͬ��

GVector2D directPos(int a, int b);//�ɸ��λ�ü���ѿ������� ��Ҫʱʵ��
double dis(const GVector2D& center, int a1, int a2);
//center���ѿ������꣩��(a1,a2) ���������꣩�ľ���
//���뿼��ƽ�ƶԳ��ԣ�������Ǹ����ĵľ��롣

double dis(int a1, int b1, int a2, int b2);//������ƽ�ƶԳ���

#define OPEN(fp, filename, mode)               \
	FILE *fp;                                  \
	fopen_s(&fp, filename, mode);              \
	if (fp == NULL)                            \
	{                                          \
		printf("open %s failed.\n", filename); \
		exit(1);                               \
	}

#define FSTREAM(fp,filename,mode)								\
	std::fstream fp=std::fstream(filename,mode);				\
	if(!fp.is_open()){											\
		std::cerr<<"Cannot open file"<<filename<<std::endl;		\
		exit(1);												\
	}

//������ռ�����ݵ��������ļ�
void output_real_matrix(const gsl_matrix * m, const char* filename);
void output_real_matrix(const gsl_matrix_complex * m, const char* filename);