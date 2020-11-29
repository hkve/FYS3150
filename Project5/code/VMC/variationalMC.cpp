#include "variationalMC.hpp"

VMC::VMC(func psi, func EL, int MCCs_, double step, double omega_, double alpha_, double beta_) {
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
	// For the interactive case, the two electrons need some sort of seperations (so <E> is not inf)
	// If this is need, I will make a setInitialGuess function. Asked on Piazza but no reply :(
	//R[0] = 0.79683;
	//R[3] = -0.79683;

	step_length = step;
	MCCs = MCCs_;

	// Set random number generator
	std::mt19937 generator(clock());
	std::uniform_real_distribution<double> s(0,1);
}

void VMC::Metropolis() {
	// Calculate trial step
	for(int i = 0; i < 6; i++) {
		R_trial[i] = R[i] + step_length * (s(generator)-0.5);
	}

	/* Calculate ratio of trial step vs old step
	   Note that trailFunctions.cpp contains the absolute square of WF's
	   So this is not needed here */
	double w = waveFunction(R_trial, omega, alpha, beta)/waveFunction(R, omega, alpha, beta);

	// Calculate if the updated R should be accepted
	// WF's depend on 0 <= r < inf, so w = [0,1]
	if(w >= s(generator)) {
		for(int i = 0; i < 6; i++) {
			R[i] = R_trial[i];
		}
		// Calculate new local energy
		Energy = localEnergy(R, omega, alpha, beta);
	}
}

void VMC::Run(string filename) {
	// Calculate initial local energy
	Energy = localEnergy(R, omega, alpha, beta);

	// For exp values
	double R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += sqrt(R12);

	// Open exp values file
	ofstream file("../../data/" + filename);
	

	// Dummy counter for keeping track of what index to write
	int writeCounter = 0;

	// Write initial exp values
	WriteExpectationValues(1, file); 
	
	// Loop over all MCCs
	for(int cycle = 2; cycle <= MCCs; cycle++) {
		
		// Preform metropolis algorithm
		Metropolis();
		
		// Calculate R12, energy and energy**2 for exp values
		R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += sqrt(R12); 

		// If this cycle should be written
		if(cycle == write[writeCounter]) {
			WriteExpectationValues(cycle, file);
			writeCounter += 1;
		}
	}

	// Close outfile to be a good boy
	file.close();
}

void VMC::Logspace(int MCCs_write) {
	// For writing (kindof) logspaced valus

	// Array holding what cycles to write
	write = new int[MCCs_write];

	// Logarithmic step length for written cycles  
	double dL = log10(MCCs)/(MCCs_write-1);

	int new_index = 2; // Start at 2 since first cycle is always written
	write[0] = new_index;
	for(int i = 1; i < MCCs_write; i++) {

		// Calculate ned index
		new_index = (int) pow(10, i*dL);

		// Since logspace is stupid for small integers, always let the next index be 1 bigger than the last
		if(write[i-1] >= new_index) {
			write[i] = write[i-1] + 1;
		}

		// If spacing is no problem, use the logspaced value
		else {
			write[i] = new_index;
		}
	}	
}

void VMC::WriteExpectationValues(int cycle, ofstream& file) {
	// Retrives exp values, normalize with cycle and write to file
	double E = ExpectationValues[0]/cycle;
	double EE = ExpectationValues[1]/cycle;
	double r12 = ExpectationValues[2]/cycle;
	file << cycle << " " << E << " " << EE << " " << r12 <<endl;
}

