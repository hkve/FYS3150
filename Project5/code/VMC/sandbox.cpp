#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {	
	int MCCs = (int) pow(10, 6);
	int stableAfter = (int) pow(10, 6);
	double step = 1.5;
	double omega = 1;
	double alpha = 1;

	VMC bar(psi_T1, EL_1, MCCs, step, omega, alpha);
	bar.Run("myoutfile.dat", "log", 1000);

	// Run without saving anything under the computations, only write at the end
	/*
	VMC foo(psi_T1, EL_1_Columb, MCCs, step, omega, alpha);
	foo.Run_NoSave();
	foo.WriteFinal("dump.dat"); // NOTE: appends
	*/
	return 0;
}