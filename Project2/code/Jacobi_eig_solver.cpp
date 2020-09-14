#include <iostream>
#include <tuple>
#include "Jacobi_eig_solver.hpp"

using namespace std;

Jacobi_eig_solver::Jacobi_eig_solver(double** A, int N) {
	m_A = A;
	m_N = N;
}


void Jacobi_eig_solver::max(double* max, int* k, int* l) {
	for(int i = 0; i < m_N; i++) {
		for(int j = 0; j < m_N; j++) {
			double aa = m_A[i][j]*m_A[i][j];
			if(aa > *max * *max) {
				*max = m_A[i][j];
				*k = i; *l = j;
			}	
		}
	}

}

Jacobi_eig_solver::~Jacobi_eig_solver() {
	for(int i = 0; i < m_N; i++) {
		delete [] m_A[i];
	}
	delete [] m_A;
}

void Jacobi_eig_solver::Print_matrix() {
	for(int i = 0; i < m_N; i++) {
		for(int j = 0; j < m_N; j++) {
			cout << m_A[i][j] << " ";
		}
		cout << endl;
	}
}