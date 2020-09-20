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
	void doJacobiRotation_(int k, int l);
	double** Solve();
	// void PrintEigenvalues();


	void PrintMatrix();
	~JacobiEigSolver();
};

#endif 