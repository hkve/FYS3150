#include <iostream>
#include <tuple>
#include <cmath>
#include "JacobiEigSolver.hpp"

using namespace std;

JacobiEigSolver::JacobiEigSolver(double** A, int N) {
	A_ = A;
	N_ = N;
	threshold_ = 1e-10;
}

void JacobiEigSolver::setA(double** A, int N) {
	delete [] A_;
	N_ = N;
	double** A_ = new double* [N];
	for (int i=0; i<N; i++) {
		A_[i] = new double [N];
		for (int j=0; j<N; j++) {
			A_[i][j] = A[i][j];
		}
	}
}

void JacobiEigSolver::CleanA(double threshold) {
	for (int i=0; i<N_; i++) {
		for (int j=0; j<N_; j++) {
			if (fabs(A_[i][j]) <= threshold) {
				A_[i][j] = 0.0;
			}
		}
	}
}

void JacobiEigSolver::getMax_(double* pmax, int* pk, int* pl) {
	// assuming a symmetrix matrix, such that we only need to search the upper half
	// of the matrix for the largest non-diagonal element.
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
	if (fabs(*pmax) < threshold_) {
		RUN = false;
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
	this->ComputeSC_(k, l, &c, &s);

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
	this->ComputeSC_(k, l, &c, &s);

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
	return B;
}


double** JacobiEigSolver::Solve() {
	bool RUN = true;
	double max;
	int k, l;

	double** B = new double* [N_];
	for (int i=0; i<N_; i++) {
		B[i] = new double [N_];
	}

	int iteration = 0;
	while (RUN) {
		cout << "solving..." << endl;
		max = 0.0;
		this->getMax_(&max, &k, &l);
		// cout << max << " " << k << " " << l << " " << endl;
		B = this->doJacobiRotation_(k, l);
		this->PrintMatrix(A_, N_);
		this->setA(B, N_);
		this->CleanA(threshold_);
		iteration += 1;
		if (iteration > 10) {
			RUN = false;
		}
		cout << endl;
	}

	return A_;
}


JacobiEigSolver::~JacobiEigSolver() {
	for(int i = 0; i < N_; i++) {
		delete [] A_[i];
	}
	delete [] A_;
}

void JacobiEigSolver::PrintMatrix(double** Matrix, int Dimension) {
	for(int i = 0; i < Dimension; i++) {
		for(int j = 0; j < Dimension; j++) {
			cout << Matrix[i][j] << " ";
		}
		cout << endl;
	}
}

