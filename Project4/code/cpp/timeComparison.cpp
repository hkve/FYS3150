#include "IsingModel2DParalell.hpp"
#include <omp.h>

int main(int argc, char* argv[]) {
	/*
	Tests timing for parallel vs series computation
	Args:
		Lstart (int) What lattice size to start at
		Lend (int) What lattice size to end at
		dL (int) Step between each lattice size
		MCCs (int) Number of MCCs
		Ntests (int) Number of tests to run for each L
		Nthreads (int) Number of threads to run (only for parallel)
		Filename (string) name for file to save
	*/
	int Lstart = atoi(argv[1]);
	int Lend = atoi(argv[2]);
	int dL = atoi(argv[3]);
	int MCCs = atoi(argv[4]);
	int Ntests = atoi(argv[5]);
	int Nthreads = atoi(argv[6]);
	string filename = argv[7];

	// As we are only testing time these dont matter
	double Tstart = 1;
	double dT = 0.1; 
	int stableMCCs = 0; 

	int NL = (Lend-Lstart)/dL;

	ofstream outfile("../data/"+filename, ios_base::app);
	for(int i = 0; i <= NL; i++) {
		#ifdef _OPENMP // If compiled with -fopenmp (runs in parallel)
		{
			int L = Lstart + i*dL;
			outfile << L << " ";
			omp_set_num_threads(Nthreads);
			double start = omp_get_wtime();
			#pragma omp parallel for 
			for(int N = 0; N < Ntests; N++) {
				double T = Tstart + N*dT; 
				IsingModel* problem = new IsingModel(L, MCCs, T, stableMCCs);
				problem->Initialize(0);
				problem->Solve();
				delete problem;	
			}
			double end = omp_get_wtime();
			double totalTime = (end-start);
			totalTime = double (totalTime/Ntests);
			outfile << totalTime << " " <<endl;
		}
		#else // If compiled without -fopenmp (runs in series)
		{
			int L = Lstart + i*dL;
			outfile << L << " ";
			clock_t start = clock();
			for(int N = 0; N < Ntests; N++) {
				double T = Tstart + N*dT;
				IsingModel* problem = new IsingModel(L, MCCs, T, stableMCCs);
				problem->Initialize(0);
				problem->Solve();
				delete problem;		
			}
			clock_t end = clock();
			double totalTime = (end-start)/(double)CLOCKS_PER_SEC; 
			totalTime = totalTime/(double)Ntests;
			outfile << totalTime << " " <<endl;
		}
		#endif
	}
	outfile.close();
}