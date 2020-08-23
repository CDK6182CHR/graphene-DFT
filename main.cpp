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
#include "integration.h"
#include<gsl/gsl_math.h>

using namespace std;

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
	cout << "EnCut=" << hbar * hbar * KCut * KCut / 2 / me / e << endl;
	cout << "N=" << N << endl;
	cout << "Nset=" << NSet << endl;
	cout << "RCount=" << RCount << endl;
	cout << "KCount=" << KCount << endl;
	cout << "E1s=" << E1s << endl;
	system_init();
	init_density();
	process_path();
	output_real_matrix(psi_1s1, "psi_1s1.dat");
	output_real_matrix(psi_1s2, "psi_1s2.dat");
	output_real_matrix(density, "density.dat");
	output_real_matrix(Vr, "Vr.dat");
	cout << "time: " << clock() - start << endl;
	return 0;
} 