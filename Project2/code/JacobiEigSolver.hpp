#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

class JacobiEigSolver
{
private:
	double** A_;
	int N_;
	void ComputeSC_(int k, int l, double* pc, double* ps);
	void getMax_(double* pmax, int* pk, int* pl);
public:
	JacobiEigSolver(double** A, int N);


	double** setSimilarityMatrix_(int k, int l);

	void PrintMatrix();
	~JacobiEigSolver();
};

#endif 