#include "problem.hpp"
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <cstdio>

using namespace std;

Problem::Problem(int p) {
	m_p = p;
	m_n = (int) pow(10, p);
	m_h = (double) 1/(m_n+1);
}

inline double Problem::f(double x) {
	return 100*exp(-10*x);
}

inline double Problem::analytical(double x) {
	return 1-(1-exp(-10)*x-exp(-10*x));
}

void Problem::Initialize() {
	m_matrix = new double*[m_n];
	m_L = new double*[m_n];
	m_U = new double*[m_n];

	m_btilde = new double[m_n];
	m_x = new double[m_n];

	double hh = m_h*m_h;

	for(int i = 0; i < m_n; i++) {
	
	m_matrix[i] = new double[m_n];
	m_x[i] = m_h*(i+1);
	m_btilde[i] = hh*f(m_x[i]);

	m_L[i] = new double[m_n];
	m_U[i] = new double[m_n];

		for(int j = 0; j < m_n; j++) {
			m_L[i][j] = 0; m_U[i][j] = 0;

			if(i == j) {
				m_matrix[i][j] = 2;
			}
			if(i-j == 1) {
				m_matrix[i][j] = -1;
			}
			if(j-i == 1) {
				m_matrix[i][j] = -1;
			}
		}
	}
}


void Problem::LU() {
	for (int i = 0; i < m_n; i++) {
		for (int k = i; k < m_n; k++) {
			double sum = 0.0;
			for (int j = 0; j < i; j++) {
				sum += (m_L[i][j]*m_U[j][k]);
			}
			m_U[i][k] = m_matrix[i][k] - sum;
		}


		for (int k = i; k < m_n; k++) {
			if (k==i) {
				m_L[i][i] = 1;
			}
			else {
				double sum = 0.0;
				for(int j = 0; j < i; j++) {
					sum += (m_L[k][j]*m_U[j][i]);
				}	
				m_L[k][i] = (m_matrix[k][i] - sum)/m_U[i][i];
			}
		}
	}	
}

void Problem::Forward_substitution() {
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
void Problem::Backward_substitution() {
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

void Problem::Write_to_file(string filename) {
	filename = filename + to_string(m_p) + ".txt";

	ifstream ifile(filename);
	if(ifile) {
		remove(filename.c_str());
	}

	ofstream outfile (filename);

	for(int i = 0; i< m_n; i++) {
		outfile << m_x[i] << "," << m_v[i] <<endl;
	}

	outfile.close();
	
}

void Problem::Print_problem() {
	for(int i = 0; i < m_n; i++) {
		for(int j = 0; j < m_n; j++) {
			cout << m_matrix[i][j] << " ";
		}
		cout  << m_btilde[i] <<endl;
	}
}

void Problem::Print_LU() {
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

void Problem::Delete() {
 	delete [] m_matrix;
 	delete [] m_x;
 	delete [] m_btilde;
 	delete [] m_L;
 	delete [] m_U;
 	delete [] m_y;
 	delete [] m_v;
}