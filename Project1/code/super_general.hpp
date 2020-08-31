#ifndef SUPER_GENERAL_HPP
#define SUPER_GENERAL_HPP

#include <string>

class S_general {
/*
	S_general class, everything is declared here
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
	double* m_v; // v vector

public:
	S_general(int p); // Constructor
	inline double f(double x); 
	inline double analytical(double x);
	void Initialize();
	void LU();
	void Forward_sub();
	void Backward_sub();
	void Write_to_file(std::string filename);

	void Print_problem();
	void Print_LU();
	void Print_sol();
	void Temp();
	void Delete();
};

#endif