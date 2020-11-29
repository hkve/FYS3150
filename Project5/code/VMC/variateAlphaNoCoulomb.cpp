#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {
	if(argc < 4) {
		cerr << "Bad usage, enter log10(MCCs), alphaStart, alphaEnd, dAlpha" <<endl;
		exit(1);
	}

	int MCCs = pow(10, atoi(argv[1]));
	double alphaStart = atof(argv[2]);
	double alphaEnd = atof(argv[3]);
	double dAlpha = atof(argv[4]);
	string filename = "varAlphaNoCoulomb";

	double omega = 1;
	double alpha = alphaStart;
	double step = 1.4;
	int N_alpha = int ((alphaEnd-alphaStart)/dAlpha) + 1;

	for(int i = 0; i < N_alpha; i++) {
		alpha = alphaStart + i*dAlpha;
		VMC* problem = new VMC(psi_T1, EL_1, MCCs, step, omega, alpha);
		problem->Logspace(1000);
		problem->Run(filename + to_string(i) + ".dat");
		cout << "Done alpha = " << alpha << endl; 
	}
}