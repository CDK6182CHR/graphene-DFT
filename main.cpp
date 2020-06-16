#include<iostream>
#include"base.h"
#include"k-const.h"
#include<functional>
#include<cstdlib>
#include<time.h>
#include<stdio.h>
#include<stdint.h>
#include "solve-secular.h"
#include "path-gen.h"
#include<gsl/gsl_errno.h>


using namespace std;

int main1()
{
	GVector2D k(0, 0);
	init_density();
	cal_k_consts(k);
	solve_k(k);
	gsl_matrix_complex* m = gsl_matrix_complex_alloc(RCount, RCount);
	for(int i=0;i<RCount;i++)
		for (int j = 0; j < RCount; j++) {
			GComplex Veff1_debug(int a, int b, double p);
			//gsl_matrix_set(m, i, j, VExtLocal(dis(r1,i,j))+VExtLocal(dis(r2,i,j)));
			gsl_matrix_complex_set(m, i, j, Veff1_debug(i,j,gsl_matrix_get(density,i,j)));
		}
	output_real_matrix(density, "density.dat");
	output_real_matrix(m, "Veff10.dat");
	output_real_matrix(phi_1s1, "phi_1s1.dat");
	output_real_matrix(Vr, "Vr.dat");
	for (int i = 0; i < NSet; i++)
		cout << kpsi_1s1[i] << endl;
	cout << "Int of psi_c" << endl;
	double x = 0;
	for (int i = 0; i < RCount; i++)
		for (int j = 0; j < RCount; j++) {
			GVector2D&& r = directPos(i, j);
			const GVector2D& Kh = B1;
			x += gsl_complex_abs(gsl_complex_exp(gsl_complex_rect(0,k*r+Kh*r)));
		}
	x *= sqrt(Omega) / RCount / RCount;
	cout << x << endl;

	cout << "Int of phi_c" << endl;
	double y = 0;
	for (int i = 0; i < RCount; i++)
		for (int j = 0; j < RCount; j++)
			y += gsl_matrix_get(phi_1s1, i, j);
	y *= Omega / RCount / RCount;
	cout << y << endl;
	output_reciprocal_matrix(H, "H.txt");
	return 0;
}

void handler(const char* reason,
	const char* file,
	int line,
	int gsl_errno)
{
	cout << "Handler" << gsl_errno << endl;
}

int main()
{
	gsl_set_error_handler(handler);
	clock_t start = clock();
	cout << "N=" << N << endl;
	cout << "Nset=" << NSet << endl;
	cout << "RCount=" << RCount << endl;
	cout << "KCount=" << KCount << endl;
	init_density();
	/*GVector2D k;
	k = (B1 * 0.2) + (B2 * 0.3);
	solve_k(k);
	for (int i = 0; i < NSet; i++) {
		cout << "E[" << i << "]=" << gsl_vector_get(Ek, i) << endl;
	}*/
	process_path();
	output_real_matrix(psi_1s1, "psi_1s1.dat");
	output_real_matrix(psi_1s2, "psi_1s2.dat");
	output_real_matrix(density, "density.dat");
	output_real_matrix(Vr, "Vr.dat");
	cout << "time: " << clock() - start << endl;
	return 0;
} 