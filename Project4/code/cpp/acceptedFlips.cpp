#include "IsingModel2D.hpp"

int main(int argc, char const *argv[]) 
{	
	double T = atof(argv[1]);
	int L = atoi(argv[2]);
	int MCS = atof(argv[3]);
	string path = "../data/";
	string filename = path + argv[4];

	IsingModel* problem = new IsingModel(L, MCS, 1, T);
	problem->Initialize(0);
	problem->Solve();
	problem->writeAcceptedFlips(filename);

	return 0;
}