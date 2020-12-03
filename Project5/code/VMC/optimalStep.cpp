#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {
	if(argc < 10) {
		cerr << "Bad usage, but thats understandle since there are a lot of parameters! Remember to include:" <<endl;
		cerr << "log(10) of MCCs" << endl;
		cerr << "omega" << endl;
		cerr << "Alpha start" <<endl; 
		cerr << "Alpha end" <<endl;
		cerr << "Alpha spacing (spacing of the alpha values to test)" << endl;
		cerr << "Step start" <<endl;  
		cerr << "Step end" <<endl;  
		cerr << "Step spacing (spacing of the steps to test)" <<endl;
		cerr << "Outfilename" <<endl;  
	}

	int MCCs = atoi(argv[1]);
	MCCs = (int)pow(10,6);
	
	double omega = atof(argv[2]); 

	double alphaStart = atof(argv[3]);
	double alphaEnd = atof(argv[4]);
	double dAlpha = atof(argv[5]);

	double stepStart = atof(argv[6]);
	double stepEnd = atof(argv[7]);
	double dStep = atof(argv[8]);

	string filename = argv[9];

	// Round since compiler somtimes truncate
	int N_alpha = (int)round((alphaEnd-alphaStart)/dAlpha) + 1; 
	int N_step = (int)round((stepEnd-stepStart)/dStep) + 1;

	double alpha; double step;
	for(int i = 0; i < N_alpha; i++) {
		alpha = alphaStart + i*dAlpha;
		for(int j = 0; j < N_step; j++) {
			step = stepStart + j*dStep;
			VMC* problem = new VMC(psi_T1, EL_1_Coulomb, MCCs, step, omega, alpha);
			problem->Run_NoSave();
			problem->WriteAccepts(filename);
		}
		cout << "Done alpha = " << alpha << endl;
	}
}