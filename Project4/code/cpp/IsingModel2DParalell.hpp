#ifndef ISINGMODEL2DPARALELL
#define ISINGMODEL2DPARALELL

#include <random> // For random numbers
#include <cstdlib> // For random device
#include <iostream> // For IO
#include <cmath>   // For exp()
#include <fstream> // For write2file
#include <omp.h> // For paralell
#include <time.h>
#include <string>
#include <iomanip>

using namespace std;

class IsingModel {
private:
	int L; // Dimension of LxL lattice
	double T; // Temperature
	int MCCs; // Number of cycles
	int stableMCCs;
	int MCCs_write; // Number of cycles to write
	int** spins; // Holding all spins
	
	double Energy, Magnetization;
	double boltzman[17]; 
	double ExpectationValues[5]; // To store E, M, E², M², |M|

	std::random_device rd;
	std::mt19937_64 generator;
	std::uniform_real_distribution<double> fdistro;
	std::uniform_int_distribution<int> idistro;

public:
	// Logic functions
	IsingModel(int L_, int MCCs_ , double T_, int stableMCCs_);
	void Initialize(int value);
	inline int PBC(int idx) {return (idx+L)%L;} // For periodic boundary conditions
	void Metropolis();
	void Solve();
	~IsingModel();

	// Calculating state variables
	int initEnergy();
	int initMagnetization();

	// Writing to file
	void writeFinalExpValues(string filename);
};

#endif