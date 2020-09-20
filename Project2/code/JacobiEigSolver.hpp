#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

class JacobiEigSolver
{
private:
	double** A_;
	double** U_;
	int N_;
	void ComputeSC_(int k, int l, double* pc, double* ps);
	bool RUN;
	double tolerance_;
	
public:
	JacobiEigSolver(double** A, int N);

	void setA(double** A, int N);
	void setTolerance(double tolerance);
	void CleanA(double tolerance);
	void getMax_(double* pmax, int* pk, int* pl);
	void doJacobiRotation_(int k, int l);
	double** Solve();
	// void PrintEigenvalues();


	void PrintMatrix(double** Matrix, int Dimension);
	~JacobiEigSolver();
};

#endif