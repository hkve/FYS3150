#ifndef JACOBI_EIG_SOLVER_HPP
#define JACOBI_EIG_SOLVER_HPP

class Jacobi_eig_solver
{
private:
	double** m_A;
	int m_N;
public:
	Jacobi_eig_solver(double** A, int N);

	void Print_matrix();
	//~Jacobi_eig_solver();
};

#endif 