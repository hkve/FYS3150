#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

class JacobiEigSolver
{
private:
	double** A_;
	double** U_;
	int N_;
	void ComputeCS_(int k, int l, double* pc, double* ps);
	bool RUN;
	double tolerance_;
	
public:
	JacobiEigSolver(double** A, int N);

	void setA(double** A, int N);
	void setTolerance(double tolerance);
	void CleanMatrix(double** A, double tolerance);
	void getMax_(double* pmax, int* pk, int* pl);
	void doJacobiRotation_(int k, int l);
	void armadilloEig();
	double** Solve();

	void PrintMatrix(double** Matrix, int Dimension);
	~JacobiEigSolver();
};

#endif