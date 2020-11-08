#ifndef ISINGMODEL2D
#define ISINGMODEL2D

#include <random> // For random numbers
#include <cstdlib> // For random device
#include <iostream> // For IO
#include <cmath>   // For exp()
#include <fstream> // For write2file

using namespace std;

class IsingModel {
private:
	int L; // Dimension of LxL lattice
	double T; // Temperature
	int MCS; // Number of cycles
	int MCS_write; // Number of cycles to write
	int** spins; // Holding all spins
	int* idx; // Index array that worked surprisingly well

	double Energy, Magnetization;
	double boltzman[17]; 
	double ExpectationValues[5]; // To store E, M, E², M², |M|

	std::mt19937_64 generator; // Random number generator 
	std::uniform_real_distribution<double> fdistro; // To initialize spins 
	std::uniform_int_distribution<int> idistro;    // To choose random spin flip
public:
	// Logic functions
	IsingModel(int L_, int MCS_ , int MCS_write_, double T_);
	void Initialize(int value);
	void Metropolis();
	void Solve();
	~IsingModel();

	// Calculating state variables
	int initEnergy();
	int initMagnetization();

	// Writing to file
	void writeLattice(ofstream& file);	

	// Functions usefull for work, migth remove
	void printSpins();

};

#endif