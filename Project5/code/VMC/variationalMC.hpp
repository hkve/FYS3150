#ifndef VARIATIONALMC_HPP
#define VARIATIONAL_HPP

#include <iostream> // Input/output
#include <cmath>	// For powers, abs
#include <string>   // For filenames
#include <fstream>  // Write to files
#include <iomanip>  // To set precision in files
#include <random>   // For uniform_real_distribution

// To not kill eyes
using namespace std; 
typedef double (*func)(double*, double, double, double);

class VMC {
private:
	// Potential and variational parameters
	double omega;
	double alpha;
	double beta;

	// Trial wavefunction and local energy
	func waveFunction;
	func localEnergy;

	// Quantom degrees of freedom and trial step
	// 6 element array, stored as [x1,y1,z1,x2,y2,z2]
	double* R; 
	double* R_trial;
	
	// Step length
	double step_length;

	// Energy and expectation values
	double Energy;
	double ExpectationValues[3];

	// Generator and distribution to calculate the ratio w >= P(R_trial)/P(R)
	std::mt19937 generator;
	std::uniform_real_distribution<double> s;

public:
	// Constructor for storing problem constants and functions
	VMC(func psi, func EL, double step, double omega_, double alpha_, double beta_=0);

	// Preform a step of the metropolis algorithm
	void Metropolis();

	// Preforms the metropolis algorithm MCCs times, how many to write and outfilename
	void Run(int MCCs, int MCCs_write, string filenmame);

	// Writing functions
	void WriteExpectationValues(int cycle, ofstream& file);
};

#endif