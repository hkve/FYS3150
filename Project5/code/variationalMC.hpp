#ifndef VARIATIONALMC_HPP
#define VARIATIONAL_HPP

#include <iostream> // Input/output
#include <cmath>	// For powers, abs
#include <fstream>  // Write to files
#include <random>   // For uniform_real_distribution

// To not kill eyes
using namespace std; 
typedef double (*func)(double*, double, double);

class VMC {
private:
	// Potential and variational parameters
	double omega;
	double alpha;

	// Trial wavefunction and local energy
	func waveFunction;
	func localEnergy;

	// Quantom degrees of freedom and trial step
	// 6 element array, stored as [x1,y1,z1,x2,y2,z2]
	double* R; 
	double* R_trial;

	// Generator and distribution to calculate the ratio w >= P(R_trial)/P(R)
	std::mt19937 generator;
	std::uniform_real_distribution<double> w;
public:
	// Constructor for storing problem constants
	VMC(double omega_, double alpha_, func psi, func EL);

	// Preform a step of the metropolis algorithm
	void Metropolis();

	// Preforms the metropolis algorithm MCCs times
	void Run(int MCCs);
	void callFunc();
};

#endif