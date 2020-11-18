#include "IsingModel2D.hpp"

int main(int argc, char const *argv[]) {
	float T = 2.4;
	int MCCs = 1000000;
	int stableAfter = 100000;
	string filename = "../data/Energy_count24.dat";

	IsingModel* problem = new IsingModel(20,MCCs,1,T);
	problem->Initialize(0);
	problem->Solve_ExpAfterStabilize(stableAfter, filename);
}