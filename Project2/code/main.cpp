#include <iostream>
#include <cmath>
#include <string>
#include "JacobiEigSolver.hpp"

#define PI 3.14159265359

using namespace std;

int main(int argc, char const *argv[])
{
	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	problem -> armadilloEig();
	problem -> Solve();

	string filename = "test";
	problem -> writeToFile(filename);

	delete problem;
	return 0;
}