#ifndef THOMAS_SINGVAL_HPP
#define THOMAS_SINGVAL_HPP

#include <string>

class Thomas_singval {
private:
	int m_p; 
	int m_n;
	double m_h;

	double* m_x;
	double* m_d;
	double* m_f;
	double* m_u;
public:
	Thomas_singval(int p);
	void Initialize();
	inline double func(double x); // u''(x)
	inline double analytical(double x); // u(x)
	void Backward_sub();

	void Write_to_file(std::string filename);
	void Print_sol();
	void Delete();
};

#endif