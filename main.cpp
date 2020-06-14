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


using namespace std;

int main()
{
	clock_t start = clock();
	cout << "N=" << N << endl;
	cout << "Nset=" << NSet << endl;
	init_density();
	/*GVector2D k;
	k = (B1 * 0.2) + (B2 * 0.3);
	solve_k(k);
	for (int i = 0; i < NSet; i++) {
		cout << "E[" << i << "]=" << gsl_vector_get(Ek, i) << endl;
	}*/
	process_path();
	output_real_matrix(density, "density.dat");
	output_real_matrix(psi_1s1, "psi_1s1.dat");
	cout << "time: " << clock() - start << endl;
} 