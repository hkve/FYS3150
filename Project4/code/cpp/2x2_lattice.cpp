#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{
	if(argc <= 3) {
		cerr << "bad usage, enter T_start, T_end and filename. (Example ./2times2.out 1 4 ../data/2x2.out";
		exit(1);
	}

	double T_start = atof(argv[1]);
	double T_end = atof(argv[2]);
	string filename = (string) argv[3];

	for(float T = T_start; T <= T_end; T+=0.1) {
		cout << T <<endl;
		IsingModel* problem = new IsingModel(2,100000,1,T);
		problem->Initialize(0);
		problem->Solve();
		problem->writeFinalExpValues(filename);
	}

	return 0;
}