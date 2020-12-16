#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char* argv[]) {
	string mode;
	if(argc < 5) {
		cerr << "Bad usage, please enter:" <<endl;
		cerr << "log(MCCs)" <<endl;
		cerr << "Omega start (lowest bound 0.01)" <<endl;
		cerr << "Omega end " <<endl;
		cerr << "N omega (number of omegas to run)" <<endl;
		cerr << "mode ('interactive' or 'noninteractive'), default: noninteractive" <<endl;
		exit(1);
	}
	else if(argc == 6) {
		mode = argv[5];
		if(mode != "interactive" && mode != "noninteractive") {
			cerr << "'"<< mode << "' is not a valid mode try 'interactive' or 'noninteractive'" <<endl;
		}
	}

	int MCCs = pow(10, atoi(argv[1]));
	double omegaStart = atof(argv[2]);
	double omegaEnd = atof(argv[3]);
	int N_omega = atoi(argv[4]);

	double dOmega = (omegaEnd-omegaStart)/(double)(N_omega-1);
	double omega;
	double alpha;
	double beta;


	string filename = "virial_" + mode + ".dat"; 

	for(int i = 0; i < N_omega; i++) {
		omega = omegaStart + i*dOmega;
		if(mode == "noninteractive") {
			alpha = 1;
			VMC* problem = new VMC(psi_T1, EL_1, MCCs, 0, omega, alpha);
			problem->Run_NoSave();
			problem->WriteVirial(filename);
			delete problem;
		}
		if(mode == "interactive") {
			alpha = 0.994;
			beta = 0.286;
			VMC* problem = new VMC(psi_T2, EL_2, MCCs, 0, omega, alpha, beta);
			problem->Run_NoSave();
			problem->WriteVirial(filename);
			delete problem;
		}

		cout << "Done omega = " << omega << " mode = " << mode <<endl;
	}
}