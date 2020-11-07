#include "IsingModel2D.hpp"

using namespace std;

IsingModel::IsingModel(int L_, int MCS_, int MCS_write_) {
	L = L_;
	MCS = MCS_;
	MCS_write = MCS_write_;

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
	// TBC, making tomorrow
}

void IsingModel::Solve() {
	// Main loop over all MCS, Solve is a bad name will figure something out
	int write = MCS/MCS_write;

	for(int cycle = 0; cycle <= MCS; cycle++) {
		
		// To write some grids for cool plots
		if(cycle % write == 0) {
			cout << "write spin gird at cycle " << cycle <<endl;
		}
	}
}

// Calculating state variables
int IsingModel::Energy() {
	int E;
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			E -= spins[i][j]*(spins[idx[i+1]][idx[j]] + spins[idx[i]][idx[j+1]]);
		}
	}
	return E;
}
int IsingModel::Magnetization() {
	int M;
	for(int i = 0; i < L; i++) {
		for(int j = 0; j < L; j++) {
			M += spins[i][j];
		}
	}
	return M;
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

int main(int argc, char const *argv[])
{
	IsingModel* problem = new IsingModel(2, 1000, 10);
	
	problem->Initialize(1);
	problem->printSpins();
	cout << "<E> = " << problem->Energy() <<endl;
	cout << "<M> = " << problem->Magnetization() <<endl;
	problem->Solve();
	delete problem;
	
	return 0;
}