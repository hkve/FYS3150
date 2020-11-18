#include "IsingModel2D.hpp"

int main(int argc, char const *argv[]) {
	float T = atof(argv[1]);
	int stableAfter = atoi(argv[2]);
	int MCCs = atoi(argv[3]);
	string filename = argv[4];

	cout << "T = " << T <<endl;
	IsingModel* problem = new IsingModel(20,MCCs,1,T);
	problem->Initialize(0);
	problem->Solve_ExpAfterStabilize(stableAfter, filename);
}