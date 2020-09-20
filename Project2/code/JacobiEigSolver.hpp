#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

class JacobiEigSolver
{
private:
	double** A_;
	int N_;
	void ComputeSC_(int k, int l, double* pc, double* ps);
	
public:
	JacobiEigSolver(double** A, int N);

	void setA(double** A);
	void getMax_(double* pmax, int* pk, int* pl);
	double** doJacobiRotation_(int k, int l);


	void PrintMatrix();
	~JacobiEigSolver();
};

#endif 