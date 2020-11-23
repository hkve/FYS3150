#include "IsingModel2D_IDX.hpp"


IsingModel::IsingModel(int L_, int MCS_, int MCS_write_, double T_) {
	L = L_;
	MCS = MCS_;
	MCS_write = MCS_write_;
	T = T_;


	// Make L+2 index array holding [L-1, 0, 2 ... L-1, 0]
	idx = new int[L+2];
	idx[0] = L-1; // Start point
	idx[L+1] = 0; // End point
	
	// Make LxL array for spins
	spins = new int*[L];
	for(int i = 0; i < L; i++) {
		spins[i] = new int[L]; // Make room for spin rows 
		idx[i+1] = i; 		   // Fill rest of idx
	}

	for(int i = 0; i < 17; i++) {boltzman[i] = 0;}
	for(int dE = -8; dE <= 8; dE += 4) {boltzman[dE+8] = exp(-dE/T);}
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;}

	std::random_device rd; 
	std::mt19937_64 generator (rd()); // Bug, different runs create same result, migth be DMA related, will investigate
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

void IsingModel::Metropolis() {
	int neighbours;
	int dE; 

	for(int i = 0; i < L*L; i++) {
		int ix = idistro(generator)%L+1;
		int iy = idistro(generator)%L+1;
		
		neighbours = spins[idx[ix]][idx[iy+1]] +
		 			 spins[idx[ix]][idx[iy-1]] +
		 			 spins[idx[ix+1]][idx[iy]] +
		 			 spins[idx[ix-1]][idx[iy]];

		
		dE = 2*spins[idx[ix]][idx[iy]]*neighbours;

		//if(dE < 0 || fdistro(generator) <= boltzman[dE + 8]) {
		if(fdistro(generator) <= boltzman[dE + 8]) {
			spins[idx[ix]][idx[iy]] *= -1;
			Energy += (double)dE;
			Magnetization += (double)2*spins[idx[ix]][idx[iy]];
		}
	}
}

void IsingModel::Solve() {
	// Main loop over all MCS, Solve is a bad name will figure something out
	int write = MCS/MCS_write;
	Energy = (double)initEnergy();
	Magnetization = (double)initMagnetization();
	
	ofstream outfile_lattice("../data/lattice.out");

	writeLattice(outfile_lattice);
	for(int cycle = 1; cycle <= MCS; cycle++) {
		Metropolis();
		
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += Magnetization;
		ExpectationValues[3] += Magnetization*Magnetization;
		ExpectationValues[4] += abs(Magnetization);

		// To write some grids for cool plots
		if(cycle % write == 0) {
			//cout << (double)cycle/MCS * 100 << "%" <<endl;
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
			E -= (double)spins[i][j]*(spins[idx[i+1]][idx[j]] + spins[idx[i]][idx[j+1]]);
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

void IsingModel::writeLattice(ofstream& file) {
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			file << spins[i][j] << " ";
		}
	}
	file << endl;
}

// Free up memory
IsingModel::~IsingModel() {
	delete [] idx;
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

int main(int argc, char const *argv[])
{

for (int L = 10; L <= 50; L+=10) { 
	for(int i = 0; i < 10; i++) {
		clock_t t;
		t = clock(); 	
		IsingModel* problem = new IsingModel(L, 100000, 1, 1);
		problem->Initialize(0); // Set random init
		problem->Solve();
		delete problem;
		double totalTime = clock() - t;
		totalTime = totalTime/(double)CLOCKS_PER_SEC;
		cout << L << " " << totalTime <<endl;
	}
}
	return 0;
}