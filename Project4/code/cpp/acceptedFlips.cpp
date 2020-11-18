#include "IsingModel2D.hpp"

int main(int argc, char const *argv[]) 
{	
	double T = atof(argv[1]);
	int L = atoi(argv[2]);
	int MCS = atof(argv[3]);
	int orientation = atoi(argv[4]);
	string filename = argv[5];

	IsingModel* problem = new IsingModel(L, MCS, 1, T);
	problem->Initialize(orientation);
	problem->Solve();
	problem->writeAcceptedFlips(filename);

	return 0;
}