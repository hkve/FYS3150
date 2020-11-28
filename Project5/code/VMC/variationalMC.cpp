#include "variationalMC.hpp"

VMC::VMC(func psi, func EL, double step, double omega_, double alpha_, double beta_) {
	// Set frequancy and variational parameter
	omega = omega_;
	alpha = alpha_;
	beta = beta_;

	waveFunction = psi;
	localEnergy = EL;

	// Create R and R_trial
	R = new double[6];
	R_trial = new double[6];

	for(int i =0; i < 6; i++) {R[i]=0; R_trial[i]=0;} // Fill R and R_trial with zeros 
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;} // Fill expectation values with zeros
	R[0] = 0.1;
	R[3] = -0.1;

	step_length = step;

	// Set random number generator
	std::mt19937 generator(clock());
	std::uniform_real_distribution<double> s(0,1);
}

void VMC::Metropolis() {
	// Implement the metropolis algorithm

	// Calculate trial step
	for(int i = 0; i < 6; i++) {
		R_trial[i] = R[i] + step_length * (s(generator)-0.5);
	}

	double w = waveFunction(R_trial, omega, alpha, beta)/waveFunction(R, omega, alpha, beta);

	if(w >= s(generator)) {
		for(int i = 0; i < 6; i++) {
			R[i] = R_trial[i];
		}

		Energy = localEnergy(R, omega, alpha, beta);
	}
}

void VMC::Run(int MCCs, int MCCs_write, string filename) {
	int write = MCCs/MCCs_write;

	Energy = localEnergy(R, omega, alpha, beta);

	double R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += sqrt(R12);

	ofstream file("../../data/" + filename);
	WriteExpectationValues(1, file);

	for(int cycle = 1; cycle <= MCCs; cycle++) {
		Metropolis();
	

		R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += sqrt(R12); 

		if(cycle%write == 0) {
			WriteExpectationValues(cycle, file);
		}
	}

	file.close();
}

void VMC::WriteExpectationValues(int cycle, ofstream& file) {
	double E = ExpectationValues[0]/cycle;
	double EE = ExpectationValues[1]/cycle;
	double r12 = ExpectationValues[2]/cycle;
	file << cycle << " " << E << " " << EE << " " << r12 <<endl;
}

