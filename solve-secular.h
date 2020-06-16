#pragma once
#include"GVector2D.h"
#include"GComplex.h"
#include<gsl/gsl_matrix.h>


extern gsl_matrix* density;  //�����ܶ�
extern gsl_matrix_complex* Vr;  //���ռ��Ƴ�
extern gsl_matrix_complex* H;  //Hamiltonian����
extern gsl_vector* Ek;  //��ǰ���ı���ֵ
extern gsl_matrix_complex* aKhs;

void init_density();  //�����ʼ������һ�Σ�������ʼ�ܶ�
bool solve_k(const GVector2D& k);  //�������k������⡣�������Ek�������Ƿ�����

void construct_H(const GVector2D& k);

void construct_Vr();

double cal_density(const GVector2D& k,int step);//��������ܶȣ�������������

double VExtLocal(double r, double p);
GComplex Veff1(int a, int b, double p); 
GComplex Veff2(int a, int b, double p);

void cal_density_single(const GVector2D& k,int n);//����һ���ܴ��϶�Ӧ���ӵ�������n: ����ֵ���

double VLDA(double p);  //LDA���������ƣ���Ϊ�ܶ�p�ķ���

GComplex V_FT(const GVector2D& Kh);  //�ֹ���ֵ���ּ���FT�����ص���V(Kh)

void solve_single_secular(const GVector2D& k);

void norm_density();// ���ܶ�����һ������