#include <iostream>
#include <cmath>
#include "JacobiEigSolver.hpp"

#define PI 3.14159265359

using namespace std;

int main(int argc, char const *argv[])
{
	int N = 3;
	double h = 1 / (N-1.0);

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

	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	
	// double max = 0;
	// int k, l;
	// problem->getMax_(&max, &k, &l);

	// cout << max << " " << k << " " << l << endl;

	// double** S = problem->setSimilarityMatrix_(k, l);
	// for(int i = 0; i < N; i++) {
	// 	for(int j = 0; j < N; j++) {
	// 		cout << S[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	problem -> Solve();

	for(int j = 1; j < N+1; j++) {
		double d = 2/(h*h); double a = -1/(h*h);
		double lambda = d + 2*a*cos(j*PI/(N+1));
		cout << "lambda = " << lambda<<endl;
		for(int i = 1; i < N+1; i++) {
			cout << sin(i*j*PI/(N+1)) << "    ";
		}
		cout <<endl<<endl;
	}

	delete problem;
	return 0;
}