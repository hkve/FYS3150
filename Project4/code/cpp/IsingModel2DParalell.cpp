#include "IsingModel2DParalell.hpp"


IsingModel::IsingModel(int L_, int MCCs_, double T_, int stableMCCs_) {
	L = L_; // Grid size is L*L
	MCCs = MCCs_; // Number of cycles to run
	stableMCCs = stableMCCs_;
	T = T_; // Temperature

	// Make LxL array for spins
	spins = new int*[L];
	for(int i = 0; i < L; i++) {
		spins[i] = new int[L]; // Make room for spin rows 
	}

	for(int i = 0; i < 17; i++) {boltzman[i] = 0;} // Fill boltzman factors to zero
	for(int dE = -8; dE <= 8; dE += 4) {boltzman[dE+8] = exp(-dE/T);} // Fill boltzman factors with actual values, indexed with dE
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;} // Fill expectation values with zeros 

	// Set random seed and generator
	std::random_device rd; 
	std::mt19937_64 generator (rd());

	// Uniform float distro for comparing with the boltzman factor
	// Uniform int distro for picking random spin 
	std::uniform_real_distribution<double> fdistro(0,1);
	std::uniform_int_distribution<int> idistro(0,L-1);
}

void IsingModel::Initialize(int value=0) {
	// Initialization of spin lattice, value chooses init state. 0 [random, defualt] * 1 [up] * -1 [down] 
	
	// If random init state
	if(value == 0) { 
		for(int i = 0; i < L; i++) {
			for(int j = 0; j < L; j++)
				spins[i][j] = (fdistro(generator) > 0.5) ? 1 : -1;
		}
	}

	// If uniform (-1[down] or 1[up]) init state
	else if(value == 1 || value == -1) { 
		for(int i = 0; i < L; i++) {
			for(int j = 0; j < L; j++)
				spins[i][j] = value;
		}
	}

	// If strange initialization
	else {
		cerr << value << " is not a valid initialization. Use 0 (random), 1 (up) or -1 (down)";
		exit(1);
	}

}


void IsingModel::Metropolis() {
	// Preforms the metropolis over the spin grid.
	int neighbours; // To hold the sum spins for all neighbours
	int dE; // To hold the change in energy

	// Loops over grid L*L times
	for(int i = 0; i < L*L; i++) {
		// Pick a random spin with index ix, iy
		int ix = idistro(generator)%L;
		int iy = idistro(generator)%L;
		
		// Calculate sum of all neighbours of spin ix, iy
		neighbours = spins[ix][PBC(iy+1)] +
					 spins[ix][PBC(iy-1)] +
					 spins[PBC(ix+1)][iy] +
					 spins[PBC(ix-1)][iy];
		
		// Calculate the change in energy if this spin should be flipped
		dE = 2*spins[ix][iy]*neighbours;

		// If the spin should be flipped
		if(fdistro(generator) <= boltzman[dE + 8]) {
			spins[ix][iy] *= -1; // Flipp it
			Energy += (double)dE;	// Update energy
			Magnetization += (double)2*spins[ix][iy]; // Update magnetizarion
		}
	}
}

void IsingModel::Solve() {
	// Main loop over all MCCs
	// Initialize energy and magnetization
	Energy = (double)initEnergy();
	Magnetization = (double)initMagnetization();
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += Magnetization;
	ExpectationValues[3] += Magnetization*Magnetization;
	ExpectationValues[4] += abs(Magnetization);


	// Loop over all cycles
	for(int cycle = 1; cycle <= MCCs; cycle++) {
		// Preform metropolis algo for L*L grid
		Metropolis();
		// Update expectation values after L*L flip attempts 
		if (cycle >= stableMCCs){
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += Magnetization;
		ExpectationValues[3] += Magnetization*Magnetization;
		ExpectationValues[4] += abs(Magnetization);
		}
	}
}

// Calculating state variables
int IsingModel::initEnergy() {
	int E = 0;
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			E -= spins[i][j]*(spins[PBC(i-1)][j] + spins[i][PBC(j-1)]);
		}
	}
	return E;
}
int IsingModel::initMagnetization() {
	int M = 0;
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			M += spins[i][j];
		}
	}
	return M;
}


void IsingModel::writeFinalExpValues(string filename) {
	// Takes a filename to write final expectation values, as well as T and MCCs

	// Normalize with the number of MCCs 
	double E = ExpectationValues[0]/(MCCs-stableMCCs);
	double E2 = ExpectationValues[1]/(MCCs-stableMCCs);
	double M = ExpectationValues[2]/(MCCs-stableMCCs);
	double M2 = ExpectationValues[3]/(MCCs-stableMCCs);
	double Mabs = ExpectationValues[4]/(MCCs-stableMCCs);
	double varE = E2-E*E;
	double varM = M2-Mabs*Mabs;

	ofstream outfile("../data/"+ filename, ios_base::app); // Appending to file
	outfile << E << " " << M << " " << E2 << " " 
			<< M2 << " " << Mabs << " " << varE << " "
			<< varM << " " << T << " " << MCCs << endl;
	outfile.close();
}


// Free up memory
IsingModel::~IsingModel() {
	for(int i = 0; i < L; i++) {
		delete [] spins[i];
	}
	delete [] spins;
}
