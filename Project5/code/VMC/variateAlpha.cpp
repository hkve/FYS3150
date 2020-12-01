#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {
	string mode;

	if(argc < 4) {
		cerr << "Bad usage, enter log10(MCCs), alphaStart, alphaEnd, dAlpha" <<endl;
		exit(1);
	}
	else if(argc == 6) {
		mode = argv[5];
		if(mode != "interactive" && mode != "noninteractive") {
			cerr << mode << " is not a valid mode, try 'interactive' or 'noninteractive'" <<endl;
			cerr << "if none is given, noninteractive is chosen" <<endl;
			exit(1);
		}
	}
	else{
		mode = "noninteractive";
	}

	int MCCs = pow(10, atoi(argv[1]));
	double alphaStart = atof(argv[2]);
	double alphaEnd = atof(argv[3]);
	double dAlpha = atof(argv[4]);
	
	string filenameMode = mode;
	filenameMode[0] = toupper(filenameMode[0]);
	string filename = "varAlpha" + filenameMode;

	func psi = psi_T1;
	func EL;
	if(mode == "noninteractive") {
		EL = EL_1;
	}
	else {
		EL = EL_1_Columb;
	}

	double omega = 1;
	double alpha = alphaStart;
	double step = 1.4;
	int N_alpha = (int)round((alphaEnd-alphaStart)/dAlpha) + 1; // Round since compiler somtimes truncate
	
	
	for(int i = 0; i < N_alpha; i++) {
		alpha = alphaStart + i*dAlpha;
		VMC* problem = new VMC(psi, EL, MCCs, step, omega, alpha);
		problem->Run(filename + to_string(i) + ".dat", "log", 1000);
		cout << "Done alpha = " << alpha << " " << mode << endl; 
	}
	
}