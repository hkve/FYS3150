#include <iostream>
#include <cmath>
#include "JacobiEigSolver.hpp"
#include <iomanip>

#include <string>
#include <fstream>

using namespace std;

JacobiEigSolver::JacobiEigSolver(double** A, int N, bool silent) {
	A_ = A; // The symmetric square matrix for which to solve
	N_ = N; // The dimension of the above square matrix
	U_ = new double* [N]; // Setting an identity matrix, to get the final eigenvectors
	for (int i = 0; i < N; i++) {
		U_[i] = new double [N];
		for (int j = 0; j < N; j++) {
			if (i == j) {
				U_[i][j] = 1.0;
			} else {
				U_[i][j] = 0.0;
			}
		}
	}
	tolerance_ = 1e-10; // Tolerance for what to interpret as 0
	if (silent) {
		LOUD = false;
	}
}

void JacobiEigSolver::setTolerance(double tolerance) {
	// Manually set tolerance.
	tolerance_ = tolerance;
}

void JacobiEigSolver::setA(double** A, int N) {
	// Manually set matrix A.
	A_ = A;
	N_ = N;
}

double** JacobiEigSolver::getA() {
	return A_;
}

void JacobiEigSolver::CleanMatrix(double** Matrix,double tolerance) {
	// Method to 'clean' matrix A by getting rid of elements smaller than the tolerance.
	for (int i=0; i<N_; i++) {
		for (int j=0; j<N_; j++) {
			if (fabs(Matrix[i][j]) <= tolerance) {
				Matrix[i][j] = 0.0;
			}
		}
	}
}

void JacobiEigSolver::getMax_(double* pmax, int* pk, int* pl) {
	// Method for finding the largest non-diagonal element of A.
	// we need only to search the upper half of the matrix for 
	// the largest non-diagonal element.
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
	}
}

void JacobiEigSolver::ComputeCS_(int k, int l, double* pc, double* ps) {
	// Method for finding the appropriate rotation coefficients of the Givens matrix
	// given the indecies k and l of the element to eliminate
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
	} else { // Failsafe if the program is asked to compute a rotation when none is needed
		*pc = 1.0;
		*ps = 0.0;
	}
}

void JacobiEigSolver::doJacobiRotation_(int k, int l) {
	// Method for performing a Jacobi rotation
	double c, s, akk, all, aik, ail, uik, uil;
	this->ComputeCS_(k, l, &c, &s);

	// Changing first the diagonal elements
	akk = A_[k][k];
	all = A_[l][l];
	A_[k][k] = c*c*akk - 2.0*c*s*A_[k][l] + s*s*all;
	A_[l][l] = s*s*akk + 2.0*c*s*A_[k][l] + c*c*all;
	// Eliminating the elements desired
	A_[k][l] = 0.0;
	A_[l][k] = 0.0;
	// Changing the appropriate non-diagonal elements
	for (int i=0; i<N_; i++) {
		if (i != k && i != l) {
			aik = A_[i][k];
			ail = A_[i][l];
			A_[i][k] = c*aik - s*ail;
			A_[k][i] = A_[i][k];
			A_[i][l] = c*ail + s*aik;
			A_[l][i] = A_[i][l];
		}
		// Adjusting the eigenvectors accordingly
		uik = U_[i][k];
		uil = U_[i][l];
		U_[i][k] = c*uik - s*uil;
		U_[i][l] = c*uil + s*uik;
	}
}

double** JacobiEigSolver::Solve() {
	// Central algorithm for the iterative solving through Jacobi rotations
	bool RUN = true;
	double max;
	int k, l;

	max = tolerance_+1.0;

	while (fabs(max)>tolerance_ && iterations_ <= N_*N_*N_) {
		max = 0.0;
		this->getMax_(&max, &k, &l);
		this->doJacobiRotation_(k, l);
		iterations_ += 1;
	}
	
	this->CleanMatrix(A_, tolerance_);
	this->CleanMatrix(U_, tolerance_);
	if (LOUD) {
		cout << "Done N = " << N_ << " number of iterations = " << iterations_ <<endl;
	}
	SOLVED = true;
	return A_;
}

void JacobiEigSolver::writeToFile(string filename) { 
	filename = "data/" + filename + ".dat";

	ofstream outfile (filename, ios_base::app); // Create file

	outfile << iterations_ << " " << N_ <<endl;

	for(int i = 0; i < N_; i++) {
		for(int j = 0; j < N_; j++) {
			if(j == 0) {
				outfile << setw(20) << setprecision(8) << A_[i][i];
			}
			outfile << setw(20) << setprecision(8) << U_[i][j];
		}
		outfile << endl;
	}
	outfile.close();
}

double* JacobiEigSolver::getEigvals() {
	double* eigvecs;
	if (SOLVED) {
	eigvecs = new double [N_];
	for (int i=0; i<N_; i++){
		eigvecs[i] = A_[i][i];
	}
	} else {
		cout << "The matrix equation has not yet been solved!" << endl;
		exit(1);
	}
	return eigvecs;
}

double** JacobiEigSolver::getEigvecs() {
	double** eigvecs;
	if (SOLVED) {
	eigvecs = new double* [N_];
	for (int i=0; i<N_; i++) {
		eigvecs[i] = new double [N_];
		for (int j=0; j<N_; j++) {
			eigvecs[i][j] = U_[i][j];
		}
	}
	this -> CleanMatrix(eigvecs, tolerance_);
	} else {
		cout << "The matrix equation as not yet been solved!" << endl;
		exit(1);
	}
	return eigvecs;
}

int JacobiEigSolver::getIterations() {
	if (SOLVED) {
		return iterations_;
	} else {
		cout << "The matrix equation as not yet been solved!" << endl;
		exit(1);
	}
}

JacobiEigSolver::~JacobiEigSolver() {
	for(int i = 0; i < N_; i++) {
		delete [] A_[i];
		delete [] U_[i];
	}
	delete [] A_;
	delete [] U_;
}

void JacobiEigSolver::PrintMatrix(double** Matrix, int Dimension) {
	for(int i = 0; i < Dimension; i++) {
		for(int j = 0; j < Dimension; j++) {
			cout << setprecision(3) << setw(10) << Matrix[i][j] << " ";
		}
		cout << endl;
	}
}

