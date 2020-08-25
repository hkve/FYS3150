#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include <string>

class Problem {
private:
	int m_p;
	int m_n;
	double m_h;
	double** m_matrix;
	double** m_L;
	double** m_U;

	double* m_x;
	double* m_btilde;
	double* m_y;
	double* m_v;

public:
	Problem(int p);
	inline double f(double x); 
	inline double analytical(double x);
	void Initialize();
	void LU();
	void Forward_substitution();
	void Backward_substitution();
	void Write_to_file(std::string filename);

	void Print_problem();
	void Print_LU();
	void Print_sol();
	void Delete();
};

#endif