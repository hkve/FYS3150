#include <iostream>
#include <tuple>
#include <cmath>
#include "JacobiEigSolver.hpp"
#include <iomanip>
// #include <armadillo>

using namespace std;
// using namespace arma;

JacobiEigSolver::JacobiEigSolver(double** A, int N) {
	A_ = A;
	N_ = N;
	threshold_ = 1e-10;
}

void JacobiEigSolver::setA(double** A, int N) {
	A_ = A;
	N_ = N;
}

void JacobiEigSolver::CleanA(double threshold) {
	for (int i=0; i<N_; i++) {
		for (int j=i+1; j<N_; j++) {
			if (fabs(A_[i][j]) <= threshold) {
				A_[i][j] = 0.0;
				A_[j][i] = 0.0;
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
		cout << "shutting down" << endl;
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

void JacobiEigSolver::doJacobiRotation_(int k, int l) {
	double s, c;
	this->ComputeSC_(k, l, &c, &s);

	double akk = A_[k][k];
	double all = A_[l][l];
	A_[k][k] = akk*c*c - 2*A_[k][l]*s*c + all*s*s;
	A_[l][l] = all*c*c + 2*A_[k][l]*s*c + akk*s*s;
	A_[k][l] = 0.0;
	A_[l][k] = 0.0;
	for (int i=0; i<N_; i++) {
		if (i != k && i != l) {
			double aik = A_[i][k];
			double ail = A_[i][l];
			A_[i][k] = c*aik - s*ail;
			A_[k][i] = A_[i][k];
			A_[i][l] = c*ail + s*aik;
			A_[l][i] = A_[i][l];
		}
	}
}


double** JacobiEigSolver::Solve() {
	bool RUN = true;
	double max;
	int k, l;

	int iteration = 0;
	this->CleanA(threshold_);
	this->PrintMatrix(A_, N_);
	max = 0.0;
	this->getMax_(&max, &k, &l);

	while (fabs(max)>threshold_ && iteration < 10) {
		cout << "solving iteration " << iteration << endl;
		// cout << max << " " << k << " " << l << " " << endl;
		this->doJacobiRotation_(k, l);
		// this->setA(B, N_);
		this->CleanA(threshold_);
		this->PrintMatrix(A_, N_);
		iteration += 1;
		// if (iteration > 10) {
		// 	RUN = false;
		// }
		cout << endl;
		// this->PrintEigenvalues()
		max = 0.0;
		this->getMax_(&max, &k, &l);

	}
	return A_;
}


JacobiEigSolver::~JacobiEigSolver() {
	for(int i = 0; i < N_; i++) {
		delete [] A_[i];
	}
	delete [] A_;
}
/*
void JacobiEigSolver::PrintEigenvalues() {
	values = arga::eig_sym(A_);
	cout << "Eigenvalues are: "
	for (int=0; i<N_; i++){
		cout << values[i] << ", "
	}
	cout << endl;
}
*/
void JacobiEigSolver::PrintMatrix(double** Matrix, int Dimension) {
	for(int i = 0; i < Dimension; i++) {
		for(int j = 0; j < Dimension; j++) {
			cout << setprecision(3) << setw(10) << Matrix[i][j] << " ";
		}
		cout << endl;
	}
}

