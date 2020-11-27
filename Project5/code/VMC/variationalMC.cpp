#include "variationalMC.hpp"
#include "trialFunctions.cpp"

VMC::VMC(double omega_, double alpha_, double step, func psi, func EL) {
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

	step_length = step;

	// Set random number generator
	std::mt19937 generator(clock());
	std::uniform_real_distribution<double> s(0,1);
}

void VMC::Metropolis() {
	// Implement the metropolis algorithm

	// Calculate trial step
	for(int i = 0; i < 6; i++) {
		R_trial[i] = R[i] + step_length * s(generator);
	}

	double w = waveFunction(R_trial, omega, alpha)/waveFunction(R, omega, alpha);

	if(w >= s(generator)) {
		for(int i = 0; i < 6; i++) {
			R[i] = R_trial[i];
		}

		Energy = localEnergy(R, omega, alpha);
	}
}

void VMC::Run(int MCCs) {

Energy = localEnergy(R, omega, alpha);

for(int cycle = 1; cycle <= MCCs; cycle++) {
	Metropolis();
}
cout << Energy << endl;
	
}


int main(int argc, char *argv[]) {
	VMC* problem = new VMC(1,1, 0.01, psi_T1, EL_1);
	problem->Run(1000000);
	return 0;
}