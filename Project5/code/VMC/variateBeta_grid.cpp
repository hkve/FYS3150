#include "variationalMC.hpp"
#include "trialFunctions.cpp"

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
		exit(1);
	}
	int MCCs = (int)pow(10,atoi(argv[1]));
	double omega = atof(argv[2]);

	double minAlpha = atof(argv[3]);
	double maxAlpha = atof(argv[4]);
	int N_alpha = atoi(argv[5]);
	double dAlpha = (maxAlpha-minAlpha)/(double)(N_alpha-1);
	double alpha;

	double minBeta = atof(argv[6]);
	double maxBeta = atof(argv[7]);
	int N_beta = atoi(argv[8]);
	double dBeta = (maxBeta-minBeta)/(double)(N_beta-1);
	double beta;

	string filename = argv[9];

	for(int i = 0; i < N_alpha; i++) {
		for(int j = 0; j < N_beta; j++) {
			
			// Calculate variational parameters
			alpha = minAlpha + i*dAlpha;
			beta = minBeta + j*dBeta;
			

			// Run sim
			VMC* problem = new VMC(psi_T2, EL_2, MCCs, 0, omega, alpha, beta);
			problem->Run_NoSave();
			problem->WriteVariational(filename);
			delete problem;

			cout << "Alpha = " << alpha << " beta = " << beta << endl;
		}
	}
}