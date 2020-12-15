#include "variationalMC.hpp"
#include "trialFunctions.cpp"
#include <omp.h>

int main(int argc, char* argv[]) {
	if(argc < 10) {
		cerr << "Bad usage, but thats understandle since there are a lot of parameters! Remember to include:" <<endl;
		cerr << "log(10) of MCCs" << endl;
		cerr << "omega" << endl;
		cerr << "Alpha min" <<endl; 
		cerr << "Alpha max" <<endl;
		cerr << "Number of alpha points" << endl;
		cerr << "Beta min" <<endl;  
		cerr << "Beta max" <<endl;  
		cerr << "Number of beta points" <<endl;
		cerr << "Outfilename" <<endl;  
		cerr << "Number of threads"<<endl;
		exit(1);
	}
	int MCCs = (int)pow(10,atoi(argv[1]));
	double omega = atof(argv[2]);

	double minAlpha = atof(argv[3]);
	double maxAlpha = atof(argv[4]);
	int N_alpha = atoi(argv[5]);
	double dAlpha = (maxAlpha-minAlpha)/(double)(N_alpha-1);

	double minBeta = atof(argv[6]);
	double maxBeta = atof(argv[7]);
	int N_beta = atoi(argv[8]);
	double dBeta = (maxBeta-minBeta)/(double)(N_beta-1);

	int Nthreads = atoi(argv[9]);
	string filename = argv[10];

	double* alphas = new double[N_beta*N_beta];
	double* betas = new double[N_beta*N_beta];


	for(int i = 0; i < N_alpha; i++) {
		for(int j = 0; j < N_beta; j++) {			
			// Calculate variational parameters
			alphas[i*N_alpha+j] = minAlpha + i*dAlpha;
			betas[i*N_alpha+j] = minBeta + j*dBeta;
		}
	}

	omp_set_num_threads(Nthreads);
	#pragma omp parallel for 
	for(int i = 0; i < N_alpha*N_beta; i++) {
			// Run sim
			VMC* problem = new VMC(psi_T2, EL_2, MCCs, 0, omega, alphas[i], betas[i]);
			problem->Run_NoSave();
			problem->WriteVariational(filename, i);
			delete problem;
			cout << omp_get_thread_num() << " done, alpha= " << alphas[i] << " beta= " << betas[i] <<endl;
	}
}