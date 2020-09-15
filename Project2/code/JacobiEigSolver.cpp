#include <iostream>
#include <tuple>
#include <cmath>
#include "JacobiEigSolver.hpp"

using namespace std;

JacobiEigSolver::JacobiEigSolver(double** A, int N) {
	A_ = A;
	N_ = N;
}

void JacobiEigSolver::setA(double** A) {
	A_ = A;
}

void JacobiEigSolver::getMax_(double* pmax, int* pk, int* pl) {
	// assuming a symmetrix matrxi, such that we only need to search the upper half
	// of the matrix for the largest non-diagonal element.
	for(int i = 0; i < N_; i++) {
		for(int j = i+1; j < N_; j++) {
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

double** JacobiEigSolver::doJacobiRotation_(int k, int l) {
	double s, c;
	JacobiEigSolver::ComputeSC_(k, l, &c, &s);

	double** B = new double* [N_];
	for (int i=0; i<N_; i++) {
		B[i] = new double [N_];
	}

	for (int i=0; i<N_; i++) {
		for (int j=0; j<N_; j++) {
			if (j==k && i!=k && i!=l) {
				B[i][j] = A_[i][k] * c - A_[i][l] * s;
			}
			if (j==l && i!=k && i!=l) {
				B[i][j] = A_[i][k] * c + A_[i][l] * s;
			}
			if (i==k && j==k) {
				B[i][j] = A_[k][k]*c*c - 2*A_[k][l]*c*s + A_[l][l]*s*s;
			}
			if (i==l && j==l) {
				B[i][j] = A_[l][l]*c*c + 2*A_[k][l]*c*s + A_[k][k]*s*s;
			}
			if (i==k && j==l) {
				B[i][j] = 0.0;
			}
			else {
				B[i][j] = A_[i][j];
			}
		}
	}

	for(int i = 0; i < N_; i++) {
		for(int j = 0; j < N_; j++) {
			cout << B[i][j] << " ";
		}
		cout << endl;
	}

	return B;
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

