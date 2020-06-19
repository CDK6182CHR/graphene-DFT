#pragma once
#include"GVector2D.h"
#include"GComplex.h"
#include<gsl/gsl_matrix.h>


extern gsl_matrix* density;  //电子密度
extern gsl_matrix_complex* Vr;  //正空间势场
extern gsl_matrix_complex* H;  //Hamiltonian矩阵
extern gsl_vector* Ek;  //当前求解的本征值
extern gsl_matrix_complex* aKhs;

extern double V0;  //平均势能

void init_density();  //程序初始化调用一次，构建初始密度
bool solve_k(const GVector2D& k);  //迭代解决k点的问题。结果存入Ek。返回是否收敛

void construct_H(const GVector2D& k);

void construct_Vr();

double Vext(const GVector2D& r0, int a, int b);
inline GComplex Veff1(int a, int b, double p); 
inline GComplex Veff2(int a, int b, double p);
double Vee(int a, int b);

double cal_density(const GVector2D& k, int step);//计算电子密度，并返回总能量
void cal_density_single(const GVector2D& k,int n);//计算一条能带上对应电子的能量。n: 特征值编号

double VLDA(double p);  //LDA交换关联势，作为密度p的泛函

GComplex V_FT(const GVector2D& Kh);  //手工数值积分计算FT。返回的是V(Kh)

void solve_single_secular(const GVector2D& k);

void norm_density();// 对密度做归一化处理

void norm_density_single(gsl_matrix* m, int electrons);