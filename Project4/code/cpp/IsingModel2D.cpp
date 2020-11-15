#include "IsingModel2D.hpp"


IsingModel::IsingModel(int L_, int MCS_, int MCS_write_, double T_) {
	L = L_;
	MCS = MCS_;
	MCS_write = MCS_write_;
	T = T_;


	acceptedFlips = 0; 	
	// Make LxL array for spins
	spins = new int*[L];
	for(int i = 0; i < L; i++) {
		spins[i] = new int[L]; // Make room for spin rows 
	}

	for(int i = 0; i < 17; i++) {boltzman[i] = 0;}
	for(int dE = -8; dE <= 8; dE += 4) {boltzman[dE+8] = exp(-dE/T);}
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;}
}

void IsingModel::Initialize(int value=0) {
	// Initialization of spin lattice, value chooses init state. 0 [random, defualt] * 1 [up] * -1 [down] 
	std::random_device rd; 
	std::mt19937_64 generator (rd()); 
	std::uniform_real_distribution<double> fdistro(0,1);
	// If random init state
	if(value == 0) { 
		for(int i = 0; i < L; i++) {
			for(int j = 0; j < L; j++)
				spins[i][j] = (fdistro(generator) > 0.5) ? 1 : -1;
		}
	}

	// If uniform (-1 or 1) init state
	else if(value == 1 || value == -1) { 
		for(int i = 0; i < L; i++) {
			for(int j = 0; j < L; j++)
				spins[i][j] = value;
		}
	}
	else {
		cerr << value << " is not a valid initialization. Use 0 (random), 1 (up) or -1 (down)";
		exit(1);
	}

}


void IsingModel::Metropolis(uniform_real_distribution<double> &fdistro, 
							uniform_int_distribution<int> &idistro,
							mt19937_64 &generator) {
	int neighbours;
	int dE; 

	for(int i = 0; i < L*L; i++) {
		int ix = idistro(generator);
		int iy = idistro(generator);
		
		neighbours = spins[ix][PBC(iy+1)] +
					 spins[ix][PBC(iy-1)] +
					 spins[PBC(ix+1)][iy] +
					 spins[PBC(ix-1)][iy];
		
		dE = 2*spins[ix][iy]*neighbours;

		if(fdistro(generator) <= boltzman[dE + 8]) {
			spins[ix][iy] *= -1;
			Energy += (double)dE;
			Magnetization += (double)2*spins[ix][iy];
			acceptedFlips++;
		}
	}
}

void IsingModel::Solve() {
	// Main loop over all MCS, Solve is a bad name will figure something out
	int write = MCS/MCS_write;
	Energy = (double)initEnergy();
	Magnetization = (double)initMagnetization();
	
	ofstream outfile_lattice("../data/lattice.out");

	std::random_device rd; 
	std::mt19937_64 generator (rd()); 
	std::uniform_real_distribution<double> fdistro(0,1);
	std::uniform_int_distribution<int> idistro(0,L-1);

	writeLattice(outfile_lattice);
	for(int cycle = 1; cycle <= MCS; cycle++) {
		Metropolis(fdistro, idistro, generator);
		
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += Magnetization;
		ExpectationValues[3] += Magnetization*Magnetization;
		ExpectationValues[4] += abs(Magnetization);

		// To write some grids for cool plots
		if(cycle % write == 0) {
			writeLattice(outfile_lattice);
		}
	}
	outfile_lattice.close();
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
	double E = ExpectationValues[0]/MCS;
	double E2 = ExpectationValues[1]/MCS;
	double M = ExpectationValues[2]/MCS;
	double M2 = ExpectationValues[3]/MCS;
	double Mabs = ExpectationValues[4]/MCS;
	double varE = E2-E*E;
	double varM = M2-M*M;
	ofstream outfile(filename, ios_base::app); // Appending to file
	outfile << E << " " << M << " " << E2 << " " 
			<< M2 << " " << Mabs << " " << varE << " "
			<< varM << " " << T << " " << MCS << endl;
	outfile.close();
}

void IsingModel::writeAcceptedFlips(string filename) {
	ofstream outfile(filename, ios_base::app); // Apennding to file
	int attemptedFlips = L*L*MCS;
	double ratioFlips = (double) acceptedFlips/attemptedFlips;
	outfile << acceptedFlips << " " << attemptedFlips << " " << ratioFlips << " " << L << " " << MCS <<endl;
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
	double temp = 1/(double)MCS;
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