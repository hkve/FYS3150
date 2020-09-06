#ifndef LU_HPP
#define LU_HPP

#include <string>

class LU {
/*
	LU class, everything is declared here
*/
private:
	int m_p; // exponent
	int m_n; // matrix dims
	double m_h; // step length
	double** m_matrix; // A matrix
	double** m_L; // L matrix
	double** m_U; // U matrix

	double* m_x; // x vector
	double* m_btilde; // b_tilde vector
	double* m_y; // y vector 
	double* m_u; // v vector

public:
	LU(int p); // Constructor
	inline double f(double x); 
	inline double analytical(double x);
	void Initialize();
	void Decomp();
	void Forward_sub();
	void Backward_sub();
	void Write_to_file(std::string filename);

	void Print_problem();
	void Print_Decomp();
	void Print_sol();
	void Temp();
	void Delete();
};

#endif