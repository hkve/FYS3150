#ifndef ISINGMODEL2D
#define ISINGMODEL2D

#include <random> // For random numbers
#include <cstdlib> // For random device
#include <iostream> // For IO

class IsingModel {
private:
	int L; // Dimension of LxL lattice
	int MCS; // Number of cycles
	int MCS_write; // Number of cycles to write
	int** spins; // Holding all spins
	int* idx; // Index array that worked surprisingly well

	float ExpectationValues[5]; // To store E, M, E², M², |M|

	std::mt19937_64 generator; // Random number generator 
	std::uniform_real_distribution<double> fdistro; // To initialize spins 
	std::uniform_int_distribution<int> idistro;    // To choose random spin flip
public:
	// Logic functions
	IsingModel(int L_, int MCS_ , int MCS_write_);
	void Initialize(int value);
	void Metropolis();
	void Solve();
	~IsingModel();

	// Calculating state variables
	int Energy();
	int Magnetization();

	// Functions usefull for work, migth remove
	void printSpins();

};

#endif