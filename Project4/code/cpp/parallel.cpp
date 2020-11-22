#include "IsingModel2DParalell.hpp"
#include <omp.h>

int main(int argc, char* argv[]) {
	int L = atoi(argv[1]); 
	int MCCs = atoi(argv[2]);
	int stableMCCs = atoi(argv[3]);
	double Tstart = atof(argv[4]);
	double Tend = atof(argv[5]);
	double dT = atof(argv[6]);
	int Nthreads = atoi(argv[7]);
	string filename = argv[8];

	int Ntemps = (int) ((Tend-Tstart)/dT);
	double start = omp_get_wtime();

	omp_set_num_threads(Nthreads);
	#pragma omp parallel for 
	for(int i = 0; i <= Ntemps+1; i++) {
		double T = Tstart + i*dT;
		IsingModel* problem = new IsingModel(L, MCCs, T, stableMCCs);
		problem->Initialize(0);
		problem->Solve();
		problem->writeFinalExpValues(filename);
		delete problem;
	}
	double end = omp_get_wtime();
	double totalTime = end-start; 
	cout << "Done L = " << L << " T=[" << Tstart << "," << Tend << "] with dT = " << dT << " time spent " << totalTime << "s" <<endl; 
}