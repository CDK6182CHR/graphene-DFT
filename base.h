#pragma once
/*基本常量等的声明*/
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex.h>
#include "GVector2D.h"
#include<iostream>
#include<fstream>

//物理常数
#define NInnerOrbit 2  //芯电子轨道数量
extern const GVector2D A1, A2, B1, B2;//正空间和倒空间基矢
extern const GVector2D r1, r2;//两个原子的位置坐标
extern const double ABohr;  //Bohr半径0.529Angstrom
extern const double hbar;  //约化Planck常数
extern const double me;//电子质量
extern const double hbar_2me;
extern const double Omega;//正空间元胞体积（面积）

extern const double A0;//正空间晶格常数

extern const double E1s, E2s;  //原子轨道本征能量
extern const double Ec_all[NInnerOrbit];  //各个原子轨道的能量本征值，按数组形式
extern const int NValence; //价电子
extern const int Z;  //中心离子核电荷数。只算价电子
extern const double e,e2k;  //元电荷


//计算参数
extern const double KCut;//平面波截断半径
extern const double prec; //收敛的相对误差判据
extern const int MaxStep; //最大迭代步数。超过就判定不收敛
extern const int LHalfCount;  //单边元胞数。考察的元胞放在正中间。
extern const int LCount;  //沿每个基矢方向的晶胞数量。
extern const int N;//晶胞数量
extern const int KCount;//1BZ高对称点路径每段折线的K点数目
extern const int RCount;//正空间划分mesh的密度。将每一条基矢等分成多少段。


extern GVector2D* Rls;  //所有要考虑的正空间格点


extern const gsl_matrix* phi_1s1, * phi_1s2, * phi_2s1, * phi_2s2; //局域芯电子波函数打表

#define KPOINTS 4  //倒空间路径结点数
extern const GVector2D KPath[KPOINTS];

#define NNeigh 8
extern const GVector2D neigh[NNeigh];  //近邻矩阵

extern const gsl_matrix
* Vext_1, * Vext_2, //外场势能是个常数
* sum_ee;           //求和常数打表


double phi_1s(double r);  //1s轨道波函数，各向同性
double phi_2s(double r);  //2s轨道波函数，各向同性

GVector2D directPos(int a, int b);//由格点位置计算笛卡尔坐标 需要时实现

double simpleDis(const GVector2D& center, int a1, int a2);//不考虑平移对称性

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

//输出正空间的数据到二进制文件
void output_real_matrix(const gsl_matrix * m, const char* filename);
void output_real_matrix(const gsl_matrix_complex * m, const char* filename);
void output_reciprocal_matrix(const gsl_matrix_complex * m, const char* filename);

inline double Vee_sum_arg(int da, int db)
{
	return gsl_matrix_get(sum_ee, da + RCount - 1, db + RCount - 1);
}