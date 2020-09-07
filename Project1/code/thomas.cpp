/*
This class numerically solves the matrix equation Au = f where
A is a tridiagonal Toeplitz matrix. The right hand side is given by
the function f(x) = 100 exp(-10x), but could easily be generalized to
any arbitrary function. The boundary conditions are set as f(0)=f(1)=0,
and x € [0,1].

We will solve the equation with LU-decomposition and using 10^p points.
The class is called with an upper limit for p, and will solve the equations
for every whole number up to, and including p.
*/

#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include "thomas.hpp"


using namespace std; // To not hurt eyes

Thomas::Thomas(int p) {
	m_pmax = p; // (not sponsored) maximum power of ten to solve for
	m_n = (int) pow(10, p);
	m_h = (double) 1/(m_n+1);

	m_u = new double[m_n]; m_x = new double[m_n]; m_f = new double[m_n]; // solution u, argument x and rhs f
	m_a = new double[m_n]; m_b = new double[m_n]; m_c = new double[m_n]; // Diagonal vectors
	m_b_tilde = new double[m_n]; m_f_tilde = new double[m_n]; // Vectors for reduced matrix
}

inline double Thomas::func(double x) { // f(x)
	return 100*exp(-10*x);
}

inline double Thomas::analytical(double x) { // analytical solution u(x)
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void Thomas::Initialize() {
	/*This initializes the class with and creates vectors a, b, c
	containing the diagonal elements of our matrix A, for our
	specific problem.*/
	double hh = m_h*m_h; // Reduce FLOPS
	
	m_a[0] = m_c[m_n-1] = 0; // Useless elements of a and c
	m_c[0] = -1;
	m_b[0] = 2; // Fixing endpoints
	m_f[0] = hh*func(m_x[0]);  // Fixing endpoints

	for(int i = 1; i < m_n; i++) {
		m_x[i-1] = m_h*i;
		m_a[i] = m_c[i] = -1; m_b[i] = 2;
		m_f[i-1] = hh*func(m_x[i-1]);
	} 
	m_x[m_n-1] = m_n*m_h;
	m_f[m_n-1] = hh*func(m_x[m_n-1]);
}

void Thomas::Forward_sub() { // Reducing a elements from array
	m_b_tilde[0] = m_b[0];
	m_f_tilde[0] = m_f[0];

	double q;

	for(int i = 1; i < m_n; i++) {
		q = m_a[i]/m_b_tilde[i-1];
		m_b_tilde[i] = m_b[i] - m_c[i-1]*q;
		m_f_tilde[i] = m_f[i] - m_f_tilde[i-1]*q;
	}
}

void Thomas::Backward_sub() { // Solving the reduced array
	m_u[m_n-1] = m_f_tilde[m_n-1]/m_b_tilde[m_n-1];
	for(int i = m_n-1; i > 0; i--) {
		m_u[i-1] = (m_f_tilde[i-1]-m_c[i-1]*m_u[i])/m_b_tilde[i-1];
	}
}

void Thomas::Write_to_file(string filename) {
	filename = "data/"+ filename + to_string(m_pmax) + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	for(int i = 0; i < 4; i++) {outfile << setw(15) << setprecision(8) << 0;} // Startpoints
	outfile << endl;
	for(int i = 0; i<m_n; i++) {
		double u_exact = analytical(m_x[i]);
		outfile << setw(15) << setprecision(8) << m_x[i];
		outfile << setw(15) << setprecision(8) << m_u[i]; 
		outfile << setw(15) << setprecision(8) << analytical(m_x[i]);
		outfile << setw(15) << setprecision(8) << fabs((u_exact-m_u[i])/u_exact) <<endl;
	}
	outfile << setw(15) << setprecision(8) << 1;
	for(int i = 0; i < 3; i++) {outfile << setw(15) << setprecision(8) << 0;} // Endpoints

	outfile.close();
}

void Thomas::Print_sol() {
	cout << "x\tu" <<endl<< "0 0" <<endl; // Making sure x(0) = 0, u(0) = 0
	for(int i = 0; i < m_n; i++) {
		cout << m_x[i] << " " << m_u[i] <<endl;
	}
	cout << "1 0" <<endl;  // Making sure x(n+1) = 1, u(n+1) = 0 
}

void Thomas::Delete() { // Free up memory
	delete [] m_a; delete [] m_b; delete [] m_c;
	delete [] m_x; delete [] m_u; delete [] m_f;
	delete [] m_b_tilde; delete [] m_f_tilde;
}