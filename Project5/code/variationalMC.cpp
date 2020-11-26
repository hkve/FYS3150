#include "variationalMC.hpp"
#include "trialFunctions.cpp"

VMC::VMC(double omega_, double alpha_, func psi, func EL) {
	// Set frequancy and variational parameter
	omega = omega_;
	alpha = alpha_;

	waveFunction = psi;
	localEnergy = EL;

	// Create R and R_trial
	R = new double[6];
	R_trial = new double[6];
	// Fill with zeros
	for(int i =0; i < 6; i++) {R[i]=0; R_trial[0]=0;} 

	// Set random number generator
	std::mt19937 generator(clock());
	std::uniform_real_distribution<double> w(0,1);
}

void VMC::Metropolis() {
	// Implement the metropolis algorithm
}

void VMC::Run(int MCCs) {

for(int cycle = 1; cycle <= MCCs; cycle++) {
	// Call metropolis and update expectation values
}

}

void VMC::callFunc() {
	waveFunction(R, omega, alpha);
	localEnergy(R, omega, alpha);
}

int main(int argc, char *argv[]) {
	VMC* problem = new VMC(1,1, psi_T1, EL_1);
	problem->callFunc();
	return 0;
}