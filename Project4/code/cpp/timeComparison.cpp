#include "IsingModel2DParalell.hpp"
#include "IsingModel2D.hpp"
#include <omp.h>

int main(int argc, char* argv[]) {
	/*
	Tests timing for parallel vs series computation
	Args:
		Lstart (int) What lattice size to start at
		Lend (int) What lattice size to end at
		dL (int) Step between each lattice size
		MCCs (int) Number of MCCs
		
	*/
	int Lstart = atoi(argv[1]);
	int Lend = atoi(argv[2]);
	int dL = atoi(argv[3]);
	int MCCs = atoi(argv[4]);
	int Ntests = atoi(argv[5]);
	string filename = "time";

	int NL = (Lend-Lstart)/dL;

	for(int i = 0; i <= NL; i++) {
		#ifdef _OPENMP 
		{
			omp_set_num_threads(4);
			#pragma omp parallel for 
			for(int N = 0; N < Ntests; N++) {
				int L = Lstart + i*dL;
				cout << L <<endl;	
			}
		}
		#else 
		{
			for(int N = 0; N < Ntests; N++) {
				int L = Lstart + i*dL;
				cout << L <<endl;	
			}
		}
		#endif
	}
}