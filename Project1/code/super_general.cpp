#include "super_general.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdio>

using namespace std;

S_general::S_general(int p) {
	// Constructor, setting pow, n_steps, and step length
	m_p = p;
	m_n = (int) pow(10, p);
	m_h = (double) 1/(m_n+1);
}

inline double S_general::f(double x) { // v''(x)
	return 100*exp(-10*x);
}

inline double S_general::analytical(double x) {
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void S_general::Initialize() {
	// Making m_n long arrays
	m_matrix = new double*[m_n]; // Matrix A
	m_L = new double*[m_n]; // Lower matrix for LU
	m_U = new double*[m_n]; // Upper matrix for LU

	m_btilde = new double[m_n]; // Result of Av
	m_x = new double[m_n]; // x points

	double hh = m_h*m_h; // Reduce flops

	for(int i = 0; i < m_n; i++) {
	
	// Each element in array is filled with other arrays
	m_matrix[i] = new double[m_n];
	m_L[i] = new double[m_n];
	m_U[i] = new double[m_n];
	
	m_x[i] = m_h*(i+1); // fill x
	m_btilde[i] = hh*f(m_x[i]); // fill b_tilde


		for(int j = 0; j < m_n; j++) {
			m_L[i][j] = 0; m_U[i][j] = 0; // Fill L and U with 0

			if(i == j) {
				m_matrix[i][j] = 2; // diagonal 
			}
			if(i-j == 1) {
				m_matrix[i][j] = -1; // lower diagonal 
			}
			if(j-i == 1) {
				m_matrix[i][j] = -1; // upper diagonal
			}
		}
	}
}


void S_general::LU() {
	// Find the LU decomp with Doolittle algorithm 
	double sum;
	for (int i = 0; i < m_n; i++) { // Looping over each row
		for (int k = i; k < m_n; k++) { // Finding U 
			sum = 0.0;
			for (int j = 0; j < i; j++) {
				sum += (m_L[i][j]*m_U[j][k]);
			}
			m_U[i][k] = m_matrix[i][k] - sum;
		}


		for (int k = i; k < m_n; k++) { // Finding L
			if (k==i) {
				m_L[i][i] = 1; // Setting diagonal of L = 1 for convenience
			}
			else {
				sum = 0.0;
				for(int j = 0; j < i; j++) {
					sum += (m_L[k][j]*m_U[j][i]);
				}	
				m_L[k][i] = (m_matrix[k][i] - sum)/m_U[i][i];
			}
		}
	}	
}

void S_general::Forward_sub() {
	// Solving Ly = b_tilde
	m_y = new double[m_n];
	m_y[0] = m_btilde[0];
	for(int i = 1; i < m_n; i++) {
		double sum = 0;
		for(int j = 0; j < i; j++) {
			sum += m_L[i][j]*m_y[j];
		}
		m_y[i] = m_btilde[i] - sum;
	}
} 
void S_general::Backward_sub() {
	//Solving Uv = y 
	m_v = new double[m_n];
	int n = m_n-1;
	m_v[n] = m_y[n]/m_U[n][n];
	for(int i = n-1; i >= 0; i--) {
		double sum = 0;
		for(int j = i+1; j < n; j++) {
			sum += m_U[i][j]*m_v[j];
		}
		m_v[i] = (m_y[i]-sum)/m_U[i][i];
	} 
}

void S_general::Write_to_file(string filename) {
	filename = "data/" + filename + to_string(m_p) + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	outfile << "0, 0, 0" <<endl; // Making sure x(0) = 0, v(0) = 0
	for(int i = 0; i< m_n; i++) {
		outfile << m_x[i] << "," << m_v[i] << "," << analytical(m_x[i]) <<endl;
	}
	outfile << "1, 0, 0" <<endl; // Making sure x(n+1) = 1, v(n+1) = 0 

	outfile.close();
}

void S_general::Print_problem() {
	// Printing the A and b_tilde in Av = b_tilde
	for(int i = 0; i < m_n; i++) {
		for(int j = 0; j < m_n; j++) {
			cout << m_matrix[i][j] << " ";
		}
		cout  << m_btilde[i] <<endl;
	}
}

void S_general::Print_LU() {
	// Printing the LU decomp
	for(int i = 0; i < m_n; i++) {
		for(int j = 0; j < m_n; j++) {
			cout << m_L[i][j] << " ";
		}
		cout << endl;
	}  
	cout <<endl;
	for(int i = 0; i < m_n; i++) {
		for(int j = 0; j < m_n; j++) {
			cout << m_U[i][j] << " ";
		}
		cout << endl;
	}  
}

void S_general::Print_sol() {
	// Printing solution (v)
	cout << "x \t v(x) for n = " << m_n <<endl;
	cout << 0 << " " << 0 <<endl; // Making sure x(0) = 0, v(0) = 0
	for(int i = 0; i < m_n; i++) {
		cout << m_x[i] << " " << m_v[i] <<endl;
	}
	cout << 1 << " " << 0 <<endl; // Making sure x(n+1) = 1, v(n+1) = 0 
}

void S_general::Delete() {
	// Free up space for next exponential
 	delete [] m_matrix;
 	delete [] m_L;
 	delete [] m_U;

 	delete [] m_x;
 	delete [] m_btilde;
 	delete [] m_y;
 	delete [] m_v;
}