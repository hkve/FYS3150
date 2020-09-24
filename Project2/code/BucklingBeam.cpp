#include <iostream>
#include <cmath>
#include <string>
#include "JacobiEigSolver.hpp"

#define PI 3.14159265359

using namespace std;

void setA(double** A, int N);

int main(int argc, char** argv)
{
	bool time;
	if(argc <= 1) {
		cout << "Bad usage, enter matrix dim NxN, example ./BucklingBeam 20";
		exit(1);
	}

	int	N = atoi(argv[1]);

	double** A = new double*[N];
	for(int i = 0; i < N; i++) {
		A[i] = new double[N];	
	}
	
	setA(A, N);

	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	problem -> Solve();
	
	string filename = "BucklingBeam";
	problem -> writeToFile(filename);
	
	return 0;
}

void setA(double** A, int N) {
	double h = 1 / (double) (N+1);
	double a = -1 / (h*h);
	double d = 2 / (h*h);
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i==j) {
				A[i][j] = d;
			} else if (abs(i-j) == 1) {
				A[i][j] = a;
			} else {
				A[i][j] = 0.0;
			}
		}
	}
}