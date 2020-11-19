#include "IsingModel2D.hpp"


IsingModel::IsingModel(int L_, int MCCs_, int MCCs_write_, double T_) {
	L = L_; // Grid size is L*L
	MCCs = MCCs_; // Number of cycles to run
	MCCs_write = MCCs_write_; // Number of cycles to write (spin grids)
	T = T_; // Temperature

	acceptedFlips = 0; 	// For counting accepted flips

	// Make LxL array for spins
	spins = new int*[L];
	for(int i = 0; i < L; i++) {
		spins[i] = new int[L]; // Make room for spin rows 
	}

	for(int i = 0; i < 17; i++) {boltzman[i] = 0;} // Fill boltzman factors to zero
	for(int dE = -8; dE <= 8; dE += 4) {boltzman[dE+8] = exp(-dE/T);} // Fill boltzman factors with actual values, indexed with dE
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;} // Fill expectation values with zeros 
}

void IsingModel::Initialize(int value=0) {
	// Initialization of spin lattice, value chooses init state. 0 [random, defualt] * 1 [up] * -1 [down] 
	// Mades random device (random seed) and generator
	std::random_device rd; 
	std::mt19937_64 generator (rd());  

	// Uniform float distro in range (0,1) to place spin up or down (only used for random state ie. 0) 
	std::uniform_real_distribution<double> fdistro(0,1);
	
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


void IsingModel::Metropolis(uniform_real_distribution<double> &fdistro, 
							uniform_int_distribution<int> &idistro,
							mt19937_64 &generator) {
	// Preforms the metropolis over the spin grid.
	int neighbours; // To hold the sum spins for all neighbours
	int dE; // To hold the change in energy

	// Loops over grid L*L times
	for(int i = 0; i < L*L; i++) {
		// Pick a random spin with index ix, iy
		int ix = idistro(generator);
		int iy = idistro(generator);
		
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
			acceptedFlips++; // Count the accepted flip
		}
	}
}

void IsingModel::UpdateExpValues() {
	// Simply updates the expectation values with the energy after a sweep of the Metropolis algo
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += Magnetization;
	ExpectationValues[3] += Magnetization*Magnetization;
	ExpectationValues[4] += abs(Magnetization);
}

void IsingModel::Solve() {
	// Main loop over all MCCs, writes the spin grid to lattice.out (will remove if we don't use)
	
	// To pick out when to write spin grid
	int write = MCCs/MCCs_write;

	// Initialize energy and magnetization
	Energy = (double)initEnergy();
	Magnetization = (double)initMagnetization();
	
	// File for grid writing (useless if we don't use it)
	ofstream outfile_lattice("../data/lattice.out");

	// Set random seed and generator
	std::random_device rd; 
	std::mt19937_64 generator (rd());

	// Uniform float distro for comparing with the boltzman factor
	// Uniform int distro for picking random spin 
	std::uniform_real_distribution<double> fdistro(0,1);
	std::uniform_int_distribution<int> idistro(0,L-1);

	// Write initial lattice
	writeLattice(outfile_lattice);

	// Loop over all cycles
	for(int cycle = 1; cycle <= MCCs; cycle++) {
		// Preform metropolis algo for L*L grid
		Metropolis(fdistro, idistro, generator);
		// Update expectation values after L*L flip attempts 
		UpdateExpValues();

		// To write some grids for cool plots
		if(cycle % write == 0) {
			writeLattice(outfile_lattice);
		}
	}
	outfile_lattice.close();
}

void IsingModel::Solve_ExpAfterStabilize(int stableAfter, string filename) {
	// Main loop over all MCCs, after the loop hits "stableAfter" int, beging writing grid energy (NOT EXPECTATION) to filename

	// Initialize energy and magnetization
	Energy = (double)initEnergy();
	Magnetization = (double)initMagnetization();

	// Random seed and generator
	std::random_device rd; 
	std::mt19937_64 generator (rd());

	// Uniform float distro for comparing with the boltzman factor
	// Uniform int distro for picking random spin 	 
	std::uniform_real_distribution<double> fdistro(0,1);
	std::uniform_int_distribution<int> idistro(0,L-1);	

	// Open file for writing grid energy
	ofstream outfile_energy("../data/"+filename);

	for(int cycle = 1; cycle <= MCCs; cycle++) {
		// Preforms metropolis algo for L*L grid
		Metropolis(fdistro, idistro, generator);

		// After L*L flip attempts, update the expectation values (Useless, just for testing)
		UpdateExpValues();

		// If we have reach a stable state
		if(cycle >= stableAfter) {
			// Write the energy this energy per spin of this grid configuration
			outfile_energy << (double)Energy/L/L << endl;
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

// Write2file functions
void IsingModel::writeLattice(ofstream& file) {
	// Write spin lattice to file
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			file << spins[i][j] << " ";
		}
	}
	file << endl;
}

void IsingModel::writeFinalExpValues(string filename) {
	// Takes a filename to write final expectation values, as well as T and MCCs

	// Normalize with the number of MCCs 
	double E = ExpectationValues[0]/MCCs;
	double E2 = ExpectationValues[1]/MCCs;
	double M = ExpectationValues[2]/MCCs;
	double M2 = ExpectationValues[3]/MCCs;
	double Mabs = ExpectationValues[4]/MCCs;
	double varE = E2-E*E;
	double varM = M2-Mabs*Mabs;

	ofstream outfile("../data/"+ filename, ios_base::app); // Appending to file
	outfile << E << " " << M << " " << E2 << " " 
			<< M2 << " " << Mabs << " " << varE << " "
			<< varM << " " << T << " " << MCCs << endl;
	outfile.close();
}

void IsingModel::writeAcceptedFlips(string filename) {
	// Write the number of accepted flips, flip attempts, prob of a flip, L and the number of MCCs
	ofstream outfile("../data/"+filename, ios_base::app); // Apennding to file
	int attemptedFlips = L*L*MCCs;
	double ratioFlips = (double) acceptedFlips/attemptedFlips;
	outfile << acceptedFlips << " " << attemptedFlips << " " << ratioFlips << " " << L << " " << MCCs <<endl;
	outfile.close();
}


// Free up memory
IsingModel::~IsingModel() {
	for(int i = 0; i < L; i++) {
		delete [] spins[i];
	}
	delete [] spins;
}

// Functions usefull for work
void IsingModel::printSpins() {
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			cout << spins[i][j] << " ";
		}
		cout << endl;
	}
}

void IsingModel::printExp() {
	double temp = 1/(double)MCCs;
	cout << "<E>" << ExpectationValues[0]*temp <<endl;
	cout << "<M>" << ExpectationValues[2]*temp <<endl;
}

/*
int main(int argc, char const *argv[])
{

	IsingModel* problem = new IsingModel(2,1000,1,1);
	problem->Initialize(0);
	problem->printSpins();
	problem->Solve();
	problem->printSpins();
	return 0;
}
*/