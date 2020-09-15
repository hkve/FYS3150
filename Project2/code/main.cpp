#include <iostream>
#include <cmath>
#include <tuple>
#include "JacobiEigSolver.hpp"

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

	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	
	double max = 0;
	int k, l;
	problem->getMax_(&max, &k, &l);

	// cout << max << " " << k << " " << l << endl;

	// double** S = problem->setSimilarityMatrix_(k, l);
	// for(int i = 0; i < N; i++) {
	// 	for(int j = 0; j < N; j++) {
	// 		cout << S[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }

	double** B;
	B = problem -> doJacobiRotation_(k, l);


	max = 0.0;
	for (int i=0; i<10; i++) {
		problem -> setA(B);
		problem -> getMax_(&max, &k, &l);
		cout << max << " " << k << " " << l << " " << endl;
		B = problem -> doJacobiRotation_(k, l);
		cout << endl;
	}
	

	delete problem;
	return 0;
}