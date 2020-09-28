#include <iostream>
#include <cmath>
#include <string>
#include <armadillo>
#include <chrono>
#include "JacobiEigSolver.hpp"

using namespace std;

void setA(double** A, int N);
double armadilloEig(double** A, int N);
void write_time(double time_jacobi, double time_arma);

int main(int argc, char** argv)
{
	int time;
	if(argc <= 1) {
		cout << "Bad usage, enter matrix dim NxN, example ./BucklingBeam 20";
		exit(1);
	} else if(argc == 3) {
		time = atoi(argv[2]);
	}

	int	N = atoi(argv[1]);

	double** A = new double*[N];
	for(int i = 0; i < N; i++) {
		A[i] = new double[N];	
	}
	
	setA(A, N);

	// Only run for time taking 
	
	double time_arma;
	if(time == 1) {
		time_arma = armadilloEig(A, N);
	}

	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	auto start_jacobi = chrono::steady_clock::now();
	problem -> Solve();
	auto end_jacobi = chrono::steady_clock::now();
	double time_jacobi = chrono::duration_cast<chrono::nanoseconds>(end_jacobi - start_jacobi).count();

	
	if (time == 1) {
		write_time(time_jacobi, time_arma);
	}
	else {
		string filename = "BucklingBeam";
		problem -> writeToFile(filename);
	}

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

double armadilloEig(double **A, int N) {
	// Takes a copy of A and stores in armadillo matrix
	// Must be run BEFORE solve (since solve changes the matrix A)
	arma::mat A_ = arma::zeros(N,N);

 	for(int i = 0; i < N; i++) { // There has to be a better way to do this
 		for(int j = 0; j < N; j++) {
 			A_(i,j) = A[i][j];
 		}
	}

	// Vector for eigenvalues and matrix for eigenvectors
	arma::vec eigval;
	arma::mat eigvec;

	auto start_arma = chrono::steady_clock::now();
	// Armadillo calculating eigenvalues and eigenvectors
	arma::eig_sym(eigval, eigvec, A_);
	auto end_arma = chrono::steady_clock::now();
	double time_arma = chrono::duration_cast<chrono::nanoseconds>(end_arma - start_arma).count();

	return time_arma;
}

void write_time(double time_jacobi, double time_arma) {
	string filename = "data/time.dat";
	ofstream outfile (filename, ios_base::app); // Create file
	outfile << time_jacobi << " " << time_arma <<endl;
	outfile.close();
}