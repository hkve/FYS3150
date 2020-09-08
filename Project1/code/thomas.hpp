#ifndef THOMAS_HPP
#define THOMAS_HPP

#include <string>

class Thomas {
private:
	int m_pmax; // Max power
	int m_n; // Dims for matrix A
	double m_h;

	// Vectors u (solution) - h^2 f (u'') - x (steps)
	double* m_u;
	double* m_b;
	double* m_x;

	// Vectors a,b,c (tridiagonal elements)  
	double* m_c;
	double* m_d;
	double* m_e;

	// Vectors , b_tilde and d_tilde
	double* m_b_tilde;
	double* m_d_tilde;
public:
	Thomas(int p); // Constructor
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