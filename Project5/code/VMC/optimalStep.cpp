#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {
	double omega = 1; 
	int MCCs = (int)pow(10,6);

	double alphaStart = 0.7;
	double alphaEnd = 1.3;
	double dAlpha = 0.1;

	double stepStart = 0.5;
	double stepEnd = 20;
	double dStep = 0.5;

	// Round since compiler somtimes truncate
	int N_alpha = (int)round((alphaEnd-alphaStart)/dAlpha) + 1; 
	int N_step = (int)round((stepEnd-stepStart)/dStep) + 1;
	string filename = "optimalStep.dat";

	double alpha; double step;
	for(int i = 0; i < N_alpha; i++) {
		alpha = alphaStart + i*dAlpha;
		for(int j = 0; j < N_step; j++) {
			step = stepStart + j*dStep;
			VMC* problem = new VMC(psi_T1, EL_1_Coulomb, MCCs, step, omega, alpha);
			problem->Run_NoSave();
			problem->WriteAccepts(filename);
		}
	}
}