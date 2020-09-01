#ifndef GENERAL_HPP
#define GENERAL_HPP

#include <string>

class General {
private:
	int pmax; // Max power
	int n; // Dims for matrix A
	double h;

	// Vectors u (solution) - f (u'') - x (steps)
	double* u;
	double* f;
	double* x;

	// Vectors a,b,c (tridiagonal elements)  
	double* a;
	double* b;
	double* c;

	// Vectors , b_tilde and f_tilde
	double* b_tilde;
	double* f_tilde;
public:
	General(int p); // Constructor
	inline double func(double x); // u''(x)
	inline double analytical(double x); // u(x)
	void Initialize(); // Fill arrays
	void Forward_sub();
	void Backward_sub();

	// Writing/printing functions
	void Print_sol();
	void Write_to_file(std::string filename);

	void Delete(); // Free up memory
};

#endif