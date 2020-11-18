#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{	
	if(argc <= 4) {
		cerr << "bad usage, enter T_start, T_end, dT, MCS and outfilename.";
		exit(1);
	}

	double T_start = atof(argv[1]);
	double T_end = atof(argv[2]);
	double dT = atof(argv[3]);
	int MCS = atoi(argv[4]);
	string filename = argv[5];
	
	if(T_start == T_end) {
		IsingModel* problem = new IsingModel(2,MCS,1,T_start);
		problem->Initialize(0);
		problem->Solve();
		problem->writeFinalExpValues(filename);		
	}
	else {

	for(float T = T_start; T <= T_end; T+=dT) {
		cout << T <<endl;
		IsingModel* problem = new IsingModel(2,MCS,1,T);
		problem->Initialize(0);
		problem->Solve();
		problem->writeFinalExpValues(filename);
		}
	}

	return 0;
}