#include <iostream>
#include <tuple>
#include <cmath>
#include "JacobiEigSolver.hpp"

using namespace std;

JacobiEigSolver::JacobiEigSolver(double** A, int N) {
	A_ = A;
	N_ = N;
}

void JacobiEigSolver::getMax_(double* pmax, int* pk, int* pl) {
	for(int i = 0; i < N_; i++) {
		for(int j = 0; j < N_; j++) {
			if (i != j) {
				double aa = A_[i][j]*A_[i][j];
				if(aa > *pmax * *pmax) {
					*pmax = A_[i][j];
					*pk = i; *pl = j;
				}	
			}
		}
	}
}

void JacobiEigSolver::ComputeSC_(int k, int l, double* pc, double* ps) {
	double tau = 0.5 * (A_[l][l] - A_[k][k]) / A_[k][l];
	double tp = -tau + sqrt(1+tau*tau);
	double tm = -tau - sqrt(1+tau*tau);
	double t = (fabs(tp) < fabs(tm)) ? tp : tm;
	*pc = sqrt(1+t*t);
	*ps = t * *pc;
}

double** JacobiEigSolver::setSimilarityMatrix_(int k, int l) {
	double c, s;
	JacobiEigSolver::ComputeSC_(k, l, &c, &s);

	double** S = new double* [N_];
	for (int i=0; i<N_; i++) {
		S[i] = new double [N_];
	}

	for (int i=0; i<N_; i++) {
		for (int j=0; j<N_; j++) {
			if (i==j) {
				if (i==k) {
					S[i][j] = c;
				}
				else {
					S[i][j] = 1;
				}
			}
			if (i==k && j==l) {
				S[i][j] = s;
			}
			if (i==l && j==k) {
				S[i][j] = -s;
			}
		}
	}
	return S;
}


JacobiEigSolver::~JacobiEigSolver() {
	for(int i = 0; i < N_; i++) {
		delete [] A_[i];
	}
	delete [] A_;
}

void JacobiEigSolver::PrintMatrix() {
	for(int i = 0; i < N_; i++) {
		for(int j = 0; j < N_; j++) {
			cout << A_[i][j] << " ";
		}
		cout << endl;
	}
}

