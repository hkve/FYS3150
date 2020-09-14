#include <iostream>
#include <cmath>
#include <tuple>
#include "Jacobi_eig_solver.hpp"

using namespace std;

int main(int argc, char const *argv[])
{
	int N = 3;
	double h = 1/(double) N;

	double** A = new double*[N];
	for(int i = 0; i < N; i++) {
		A[i] = new double[N];	
	}

	double a = -1 / (h*h);
	double d = 2 / (h*h);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i==j) {
				A[i][j] = d;
			}
			if (abs(i-j) == 1) {
				A[i][j] = a;
			}
		}
	}

	Jacobi_eig_solver* problem = new Jacobi_eig_solver(A, N);

	
	double max = 0;
	int k, l;
	problem->max(&max, &k, &l);
	
	cout << max << " " << k << " " << l;

	delete problem;
	return 0;
}