#ifndef THOMAS_HPP
#define THOMAS_HPP

#include <string>

class Thomas {
private:
	int m_pmax; // Max power
	int m_n; // Dims for matrix A
	double m_h;

	// Vectors u (solution) - f (u'') - x (steps)
	double* m_u;
	double* m_f;
	double* m_x;

	// Vectors a,b,c (tridiagonal elements)  
	double* m_a;
	double* m_b;
	double* m_c;

	// Vectors , b_tilde and f_tilde
	double* m_b_tilde;
	double* m_f_tilde;
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