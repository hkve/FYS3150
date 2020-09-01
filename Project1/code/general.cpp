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
#include "general.hpp"


using namespace std; // To not hurt eyes

General::General(int p) {
	pmax = p; // (not sponsored) maximum power of ten to solve for
	n = (int) pow(10, p);
	h = (double) 1/(n+1);

	u = new double[n]; x = new double[n]; f = new double[n]; // solution u, argument x and rhs f
	a = new double[n]; b = new double[n]; c = new double[n]; // Diagonal vectors
	b_tilde = new double[n]; f_tilde = new double[n]; // Vectors for reduced matrix
}

inline double General::func(double x) { // f(x)
	return 100*exp(-10*x);
}

inline double General::analytical(double x) { // analytical solution u(x)
	return 1-(1-exp(-10))*x-exp(-10*x);
}

void General::Initialize() {
	/*This initializes the class with and creates vectors a, b, c
	containing the diagonal elements of our matrix A, for our
	specific problem.*/
	double hh = h*h; // Reduce FLOPS
	
	a[0] = c[n-1] = 0; // Useless elements of a and c
	c[0] = -1;
	b[0] = 2; // Fixing endpoints
	f[0] = hh*func(x[0]);  // Fixing endpoints

	for(int i = 1; i < n; i++) {
		x[i-1] = h*i;
		a[i] = c[i] = -1; b[i] = 2;
		f[i-1] = hh*func(x[i-1]);
	} 
	x[n-1] = n*h;
	f[n-1] = hh*func(x[n-1]);
}

void General::Forward_sub() { // Reducing a elements from array
	b_tilde[0] = b[0];
	f_tilde[0] = f[0];

	for(int i = 1; i < n; i++) {
		b_tilde[i] = b[i] - (a[i]*c[i-1])/b_tilde[i-1];
		f_tilde[i] = f[i] - (a[i]*f_tilde[i-1])/b_tilde[i-1];
	}
}

void General::Backward_sub() { // Solving the reduced array
	u[n-1] = f_tilde[n-1]/b_tilde[n-1];
	for(int i = n-1; i > 0; i--) {
		u[i-1] = (f_tilde[i-1]-c[i-1]*u[i])/b_tilde[i-1];
	}
}

void General::Write_to_file(string filename) {
	filename = "data/"+ filename + to_string(pmax) + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	outfile << "0, 0, 0" << endl; // Making sure x(0) = 0, u(0) = 0
	for(int i = 0; i< n; i++) {
		outfile << x[i] << "," << u[i] << "," << analytical(x[i]) <<endl;
	}
	outfile << "1, 0, 0" <<endl; // Making sure x(n+1) = 1, u(n+1) = 0 

	outfile.close();
}

void General::Print_sol() {
	cout << "x\tu" <<endl<< "0 0" <<endl; // Making sure x(0) = 0, u(0) = 0
	for(int i = 0; i < n; i++) {
		cout << x[i] << " " << u[i] <<endl;
	}
	cout << "1 0" <<endl;  // Making sure x(n+1) = 1, u(n+1) = 0 
}

void General::Delete() { // Free up memory
	delete [] a; delete [] b; delete [] c;
	delete [] x; delete [] u; delete [] f;
	delete [] b_tilde; delete [] f_tilde;
}
