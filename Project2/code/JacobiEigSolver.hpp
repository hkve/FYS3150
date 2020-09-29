#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

#include <string>

class JacobiEigSolver
{
private:
	bool SOLVED = false;
	bool LOUD = true;
	double** A_;
	double** U_;
	int N_;
	void ComputeCS_(int k, int l, double* pc, double* ps);
	bool RUN;
	double tolerance_;
	int iterations_ = 0;
	
public:
	JacobiEigSolver(double** A, int N, bool silent=false);

	void setA(double** A, int N);
	void setTolerance(double tolerance);
	void CleanMatrix(double** A, double tolerance);
	void getMax_(double* pmax, int* pk, int* pl);
	void doJacobiRotation_(int k, int l);
	
	void writeToFile(std::string filename);
	double** Solve();
	double** getA();
	double* getEigvals();
	double ** getEigvecs();
	int getIterations();

	void PrintMatrix(double** Matrix, int Dimension);
	~JacobiEigSolver();
};

#endif