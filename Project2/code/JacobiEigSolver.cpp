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
	U_ = new double* [N];
	for (int i=0; i<N; i++) {
		U_[i] = new double [N];
		U_[i][i] = 1.0;
	}
	// for (int i=0; i<N, i++) {
	// 	for (int j=i; j<N; j++) {
	// 		if (j==i) {
	// 			U_[i][j] = 1.0;
	// 		} else {
	// 			U_[i][j] = 0.0;
	// 		}
	// 	}
	// }
	tolerance_ = 1e-10;
}

void JacobiEigSolver::setTolerance(double tolerance) {
	tolerance_ = tolerance;
}

void JacobiEigSolver::setA(double** A, int N) {
	A_ = A;
	N_ = N;
}

void JacobiEigSolver::CleanA(double tolerance) {
	for (int i=0; i<N_; i++) {
		for (int j=i+1; j<N_; j++) {
			if (fabs(A_[i][j]) <= tolerance) {
				A_[i][j] = 0.0;
				A_[j][i] = 0.0;
			}
		}
	}
}

void JacobiEigSolver::getMax_(double* pmax, int* pk, int* pl) {
	// assuming a symmetrix matrix, such that we only need to search the upper half
	// of the matrix for the largest non-diagonal element.
	double aij;
	for(int i = 0; i < N_-1; i++) {
		for(int j = i+1; j < N_; j++) {
			aij = fabs(A_[i][j]);
			if(aij > fabs(*pmax)) {
				*pmax = A_[i][j];
				*pk = i; *pl = j;
			}	
		}
	}
	if (fabs(*pmax) < tolerance_) {
		RUN = false;
		cout <<"shutting down, max = " << fabs(*pmax) << endl;
	}
}

void JacobiEigSolver::ComputeSC_(int k, int l, double* pc, double* ps) {
	double tau, t;

	if (A_[k][l] != 0) {
		tau = 0.5 * (A_[l][l] - A_[k][k]) / A_[k][l];
		if (tau >= 0) {
			t = 1.0 / (tau + sqrt(1.0 + tau*tau));
		} else {
			t = -1.0 / (-tau + sqrt(1.0 + tau*tau));
		}
		*pc = 1.0 / sqrt(1.0 + t*t);
		*ps = *pc * t;
	} else {
		*pc = 1.0;
		*ps = 0.0;
	}
}

void JacobiEigSolver::doJacobiRotation_(int k, int l) {
	double s, c, akk, all, aik, ail, uik, uil;
	this->ComputeSC_(k, l, &c, &s);

	akk = A_[k][k];
	all = A_[l][l];
	A_[k][k] = c*c*akk - 2.0*c*s*A_[k][l] + s*s*all;
	A_[l][l] = c*c*akk + 2.0*c*s*A_[k][l] + s*s*all;
	A_[k][l] = 0.0;
	A_[l][k] = 0.0;
	for (int i=0; i<N_; i++) {
		if (i != k && i != l) {
			aik = A_[i][k];
			ail = A_[i][l];
			A_[i][k] = c*aik - s*ail;
			A_[k][i] = A_[i][k];
			A_[i][l] = c*ail + s*aik;
			A_[l][i] = A_[i][l];


		}
		uik = U_[i][k];
		uil = U_[i][l];
		U_[i][k] = c*uik - s*uil;
		U_[i][l] = c*uil + s*uik;
	}

}


double** JacobiEigSolver::Solve() {
	bool RUN = true;
	double max;
	int k, l;

	int iteration = 0;
	this->CleanA(tolerance_);
	this->PrintMatrix(A_, N_);
	max = 0.0;
	this->getMax_(&max, &k, &l);

	while (fabs(max)>tolerance_ && iteration < 100) {
		iteration += 1;
		cout << "solving iteration " << iteration << endl;
		cout << max << " " << k << " " << l << " " << endl;
		this->doJacobiRotation_(k, l);
		// this->CleanA(tolerance_);
		this->PrintMatrix(A_, N_);
		cout << endl;
		// this->PrintEigenvalues()
		max = 0.0;
		this->getMax_(&max, &k, &l);
	}
	this->PrintMatrix(U_, N_);
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

