#include "IsingModel2DParalell.hpp"
#include <omp.h>
// int main(int argc, char const *argv[])
// {   
//     // tror kommentarene under er feil, evt fjern de senere
    
//     /*args:
//     T (double): Temperature in kT/J
//     MCCs (int): Number of cycles to run 
//     initSpin (int): 0 or 1. The orientation of the init spins. 0 for random, 1 for all positively aligned
//     filename (string) Name of file to save in data folder
//     */
    
//     IsingModel* problem = new IsingModel(L,MCCs,1,T);
//     problem->Initialize(0);
//     problem->Solve(); // This saves ExpValues at every cycle
//     cout << L << " " << T << " "<< MCCs << " " << filename << endl;
//     problem->writeLattice(filename);
//     //problem->Solve();
//     //problem->writeFinalExpValues(filename);	

	

// 	return 0;
// }


int main(int argc, char* argv[]) {
	int L = atoi(argv[1]);
    double T = atof(argv[2]);
    int MCCs = atoi(argv[3]);
    int stableMCCs = atoi(argv[4]);
    string filename = (string) argv[5];

	
	double start = omp_get_wtime();

	omp_set_num_threads(4);
	//#pragma omp parallel for 
	
    IsingModel* problem = new IsingModel(L, MCCs, T, stableMCCs);
    problem->Initialize(0);
    problem->Solve();
    problem->writeLattice(filename);
    delete problem;

	double end = omp_get_wtime();
	double totalTime = end-start; 
	cout << "Done T = " << T << ", L = " << L  << ", time spent:"   << totalTime << "s" <<endl; 
    return 0;
}