#ifndef ISINGMODEL2D
#define ISINGMODEL2D

#include <random> // For random numbers
#include <cstdlib> // For random device
#include <iostream> // For IO
#include <cmath>   // For exp()
#include <fstream> // For write2file
#include <time.h>
#include <string>

using namespace std;

class IsingModel {
private:
	int L; // Dimension of LxL lattice
	double T; // Temperature
	int MCS; // Number of cycles
	int MCS_write; // Number of cycles to write
	int** spins; // Holding all spins
	
	double Energy, Magnetization;
	double boltzman[17]; 
	double ExpectationValues[5]; // To store E, M, E², M², |M|

	int acceptedFlips;
public:
	// Logic functions
	IsingModel(int L_, int MCS_ , int MCS_write_, double T_);
	void Initialize(int value);
	inline int PBC(int idx) {return (idx+L)%L;} // For periodic boundary conditions
	void Metropolis(uniform_real_distribution<double> &fdistro,
				    uniform_int_distribution<int> &idistro,  
				    mt19937_64 &generator);
	void Solve();
	~IsingModel();

	// Calculating state variables
	int initEnergy();
	int initMagnetization();

	// Writing to file
	void writeLattice(ofstream& file);	
	void writeFinalExpValues(string filename);
	void writeAcceptedFlips(string filename);

	// Functions usefull for work, migth remove
	void printSpins();
	void printExp();
};

#endif