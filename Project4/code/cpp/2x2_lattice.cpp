#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{
	if(argc <= 2) {
		cerr << "bad usage, enter T_start and T_end. (Example ./2times2.out 1 4";
		exit(1);
	}

	double T_start = atof(argv[1]);
	double T_end = atof(argv[2]);
	string filename = "../data/expValues.out";

	for(float T = T_start; T <= T_end; T+=0.1) {
		cout << T <<endl;
		IsingModel* problem = new IsingModel(2,1000,1,T);
		problem->Initialize(0);
		problem->Solve();
		problem->writeFinalExpValues(filename);
	}

	return 0;
}