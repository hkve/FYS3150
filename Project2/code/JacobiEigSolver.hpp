#ifndef JACOBIEIGSOLVER_HPP
#define JACOBIEIGSOLVER_HPP

class JacobiEigSolver
{
private:
	double** A_;
	int N_;
	void ComputeSC_(int k, int l, double* pc, double* ps);
<<<<<<< HEAD
	void getMax_(double* pmax, int* pk, int* pl);
public:
	JacobiEigSolver(double** A, int N);


=======
	
public:
	JacobiEigSolver(double** A, int N);

	void setA(double** A);
	void getMax_(double* pmax, int* pk, int* pl);
>>>>>>> f10fb064ad49117dd131a86260332c21d4b7914e
	double** setSimilarityMatrix_(int k, int l);
	double** doJacobiRotation_(int k, int l);


	void PrintMatrix();
	~JacobiEigSolver();
};

#endif 