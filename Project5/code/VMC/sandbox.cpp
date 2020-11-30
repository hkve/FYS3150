#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {	
	int MCCs = (int) pow(10, 6);
	int stableAfter = (int) pow(10, 6);
	double step = 1.5;
	double omega = 1;
	double alpha = 1.4;

	VMC noStable(psi_T1, EL_1, MCCs, step, omega, alpha);
	noStable.Run("noStable.dat");

	VMC Stable(psi_T1, EL_1, MCCs, step, omega, alpha);
	Stable.Stabilize(stableAfter);
	Stable.Run("Stable.dat");

	return 0;
}