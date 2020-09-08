#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "thomas_singval.hpp"

using namespace std;

Thomas_singval::Thomas_singval(int p) {
	m_p = p;
	m_n = (int) pow(10, p);
	
	m_h = (double) 1/(m_n+1);

	m_x = new double[m_n];
	m_d = new double[m_n];
	m_b = new double[m_n];
	m_u = new double[m_n];
}

inline double Thomas_singval::func(double x) { // u''(x)
	return 100*exp(-10*x);
}

inline double Thomas_singval::analytical(double x) {
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void Thomas_singval::Initialize() {
	double hh = m_h*m_h;
	for(int i = 0; i < m_n; i++) {
		m_x[i] = m_h*(i+1);
		m_b[i] = hh*func(m_x[i]);
		m_d[i] = (i+2)/((double) i+1);
	}
}

void Thomas_singval::Backward_sub() {
	for(int i = 1; i < m_n; i++) {
		m_b[i] = m_b[i] + m_b[i-1]/m_d[i-1]; 
	}
	
	m_u[m_n-1] = m_b[m_n-1]/m_d[m_n-1];
	for(int i = m_n-2; i >= 0; i--) {
		m_u[i] = (m_b[i]+m_u[i+1])/m_d[i];
	}
}

void Thomas_singval::Print_sol() {
	for(int i = 0; i < m_n; i++) {
		cout << m_x[i] << " " << m_u[i] <<endl;
	}
}

void Thomas_singval::Write_to_file(string filename) {
	filename = "data/"+ filename + to_string(m_p) + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	// m_m 		m_u 	analytical		rel error
	for(int i = 0; i < 4; i++) {outfile << setw(15) << setprecision(8) << 0;} // Startpoints
	outfile << endl;
	for(int i = 0; i<m_n; i++) {
		double u_exact = analytical(m_x[i]);
		outfile << setw(15) << setprecision(8) << m_x[i];
		outfile << setw(15) << setprecision(8) << m_u[i]; 
		outfile << setw(15) << setprecision(8) << u_exact;
		outfile << setw(15) << setprecision(8) << fabs((u_exact-m_u[i])/u_exact) <<endl;
	}
	outfile << setw(15) << setprecision(8) << 1;
	for(int i = 0; i < 3; i++) {outfile << setw(15) << setprecision(8) << 0;} // Endpoints


	outfile.close();
}

void Thomas_singval::Delete() { // Free up memory
	delete [] m_x;
	delete [] m_d;
	delete [] m_b;
	delete [] m_u;
}