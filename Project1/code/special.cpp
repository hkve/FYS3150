#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "special.hpp"

using namespace std;

Special::Special(int p) {
	m_p = p;
	m_n = (int) pow(10, p);
	
	m_h = (double) 1/(m_n+1);

	m_x = new double[m_n];
	m_d = new double[m_n];
	m_f = new double[m_n];
	m_u = new double[m_n];
}

inline double Special::func(double x) { // u''(x)
	return 100*exp(-10*x);
}

inline double Special::analytical(double x) {
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void Special::Initialize() {
	double hh = m_h*m_h;
	for(int i = 0; i < m_n; i++) {
		m_x[i] = m_h*(i+1);
		m_f[i] = hh*func(m_x[i]);
		m_d[i] = (i+2)/((double) i+1);
	}

	for(int i = 1; i < m_n; i++) {
		m_f[i] = m_f[i] + m_f[i-1]/m_d[i-1]; 
	}
}

void Special::Backward_sub() {
	m_u[m_n-1] = m_f[m_n-1]/m_d[m_n-1];
	for(int i = m_n-2; i >= 0; i--) {
		m_u[i] = (m_f[i]+m_u[i+1])/m_d[i];
	}
}

void Special::Print_sol() {
	for(int i = 0; i < m_n; i++) {
		cout << m_x[i] << " " << m_u[i] <<endl;
	}
}

void Special::Write_to_file(string filename) {
	filename = "data/"+ filename + to_string(m_p) + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	outfile << "0, 0, 0" << endl; // Making sure x(0) = 0, u(0) = 0
	for(int i = 0; i < m_n; i++) {
		outfile << m_x[i] << "," << m_u[i] << "," << analytical(m_x[i]) <<endl;
	}
	outfile << "1, 0, 0" <<endl; // Making sure x(n+1) = 1, u(n+1) = 0 

	outfile.close();
}

void Special::Delete() { // Free up memory
	delete [] m_x;
	delete [] m_d;
	delete [] m_f;
}