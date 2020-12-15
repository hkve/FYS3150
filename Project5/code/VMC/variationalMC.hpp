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
	
	// Step length, MCCs
	double step_length;
	int MCCs; 

	// Energy and expectation values
	double Energy;
	double ExpectationValues[5];

	// Counting the number of accepted steps
	int accepted;

	// Storing what indexes should be written to file
	int* write;

	// Generator and distribution to calculate the ratio w >= P(R_trial)/P(R)
	std::mt19937 generator;
	std::uniform_real_distribution<double> s;

public:
	// Constructor for storing problem constants and functions
	VMC(func psi, func EL, int MCCs_, double step, double omega_, double alpha_, double beta_=0);

	// Preform a step of the metropolis algorithm
	void Metropolis();

	// Constr
	void Stabilize(int stableAfter);

	// Preforms the metropolis algorithm MCCs times, how many to write and outfilename
	void Run(string filenmame, string spaced="log", int MCCs_write=1000);

	// Same as above, but nothing is saved under calculations
	void Run_NoSave();

	// Calculates the potential energy
	double potEnergy(double* r);

	// Getter for energy, used when variating both alpha and beta
	double getEnergy();

	// If step is init as 0, the optimal step is chosen
	double optimalStep();

	// Sets what mode you want to write 
	void setOutfileParameters(int MCCs_write, string spaced);

	// Different modes of writing data
	void Logspace(int MCCs_write);
	void Linspace(int MCCs_write);
	void Final();

	// Writing functions
	void WriteExpectationValues(int cycle, ofstream& file);
	void WriteVariational(string filename, int idx=-1);
	void WriteAccepts(string filename);
	void WriteVirial(string filename);
};

#endif