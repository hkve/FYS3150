#include "variationalMC.hpp"

VMC::VMC(func psi, func EL, int MCCs_, double step, double omega_, double alpha_, double beta_) {
	// Set MCCs
	MCCs = MCCs_;

	// Set frequancy and variational parameter
	omega = omega_;
	alpha = alpha_;
	beta = beta_;

	// Set class functions equal to input functions
	waveFunction = psi;
	localEnergy = EL;

	if(step == 0) {
		step_length = optimalStep();
	}
	else {
		step_length = step;
	}

	// Set random number generator
	std::mt19937 generator(clock());
	std::uniform_real_distribution<double> s(0,1);
	
	// Create R and R_trial
	R = new double[6];
	R_trial = new double[6];

	// Fill R with random step and R_trial with zeros
	for(int i =0; i < 6; i++) {
		R[i] = step_length * (s(generator)-0.5);
		R_trial[i] = 0;
	}  
	for(int i = 0; i < 5; i++) {ExpectationValues[i] = 0;} // Fill expectation values with zeros
	
	accepted = 0;
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

		// Count the accepted step
		accepted++;
	}
}

void VMC::Stabilize(int stableAfter) {
	for(int i = 0; i < stableAfter-1; i++) {
		Metropolis();
	}
}

void VMC::Run(string filename, string spaced, int MCCs_write) {
	// Choose how to space the values written to outfile
	setOutfileParameters(MCCs_write, spaced);

	// Calculate initial local energy
	Energy = localEnergy(R, omega, alpha, beta);

	// For exp values
	double R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
	R12 = sqrt(R12);
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += R12;
	ExpectationValues[3] += potEnergy(R);
	ExpectationValues[4] += 1/R12;

	// Open exp values file
	ofstream file("../data/" + filename);
	

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
		R12 = sqrt(R12);
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += R12; 
		ExpectationValues[3] += potEnergy(R);
		ExpectationValues[4] += 1/R12;

		// If this cycle should be written
		if(cycle == write[writeCounter]) {
			WriteExpectationValues(cycle, file);
			writeCounter += 1;
		}
	}

	// Close outfile to be a good boy
	file.close();
}

void VMC::Run_NoSave() {
	// Calculate initial local energy
	Energy = localEnergy(R, omega, alpha, beta);

	// For exp values
	double R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
	R12 = sqrt(R12);
	ExpectationValues[0] += Energy;
	ExpectationValues[1] += Energy*Energy;
	ExpectationValues[2] += R12;
	ExpectationValues[3] += potEnergy(R);
	ExpectationValues[4] += 1/R12;

		// Loop over all MCCs
	for(int cycle = 1; cycle <= MCCs; cycle++) {
		
		// Preform metropolis algorithm
		Metropolis();
		
		// Calculate R12, energy and energy**2 for exp values
		R12 = (R[0]-R[3])*(R[0]-R[3]) +
			     	 (R[1]-R[4])*(R[1]-R[4]) +
				     (R[2]-R[5])*(R[2]-R[5]);
		R12 = sqrt(R12);
		ExpectationValues[0] += Energy;
		ExpectationValues[1] += Energy*Energy;
		ExpectationValues[2] += R12; 
		ExpectationValues[3] += potEnergy(R);
		ExpectationValues[4] += 1/R12;
	}
}

double VMC::getEnergy() {
	// Getter for the energy
	return ExpectationValues[0]/MCCs;
}

double VMC::potEnergy(double* r) {
	double r1_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double r2_squared = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];

	return 0.5*omega*omega*(r1_squared+r2_squared);
}

double VMC::optimalStep() {
	// If step is set to zero, simply try to find the optimal step
	double m = 1.381 * pow(omega, -0.5);
	return m * pow(alpha, -0.5);
}

void VMC::setOutfileParameters(int MCCs_write, string spaced) {
	// Choosing how to space the values written to outfile
	if(spaced=="log" || spaced=="Log") {
		Logspace(MCCs_write);
	}
	else if (spaced=="lin" || spaced=="Lin") {
		Linspace(MCCs_write);
	}
	else if(spaced=="final" || spaced=="Final") {
		Final();
	}
	else {
		cerr << spaced << " is not a recognized outfile parameter. Use 'log', 'lin' or 'final'" <<endl;
		exit(1);
	}
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

void VMC::Linspace(int MCCs_write) {
	// For writing linspaced values
	
	// Linspaced index step
	int dN = MCCs/MCCs_write;

	// Array holding what cycles to write
	write = new int[MCCs_write];

	for(int i = 0; i < MCCs_write; i++) {
		write[i] = (i+1)*dN;
	}
}

void VMC::Final() {
	// For just writing the last value
	write = new int[1];
	write[0] = MCCs;
}

void VMC::WriteExpectationValues(int cycle, ofstream& file) {
	// Retrives exp values, normalize with cycle and write to file
	double E = ExpectationValues[0]/cycle;
	double EE = ExpectationValues[1]/cycle;
	double r12 = ExpectationValues[2]/cycle;
	file << cycle << " " << E << " " << EE << " " << r12 <<endl;
}

void VMC::WriteVariational(string filename) {
	// Appending parameters and expectation values to file, assuming Run or Run_NoSaved has been called
	ofstream file;
	file.open("../data/" + filename, ios_base::app);
	double E = ExpectationValues[0]/MCCs;
	double EE = ExpectationValues[1]/MCCs;
	double r12 = ExpectationValues[2]/MCCs;

	file << MCCs << " " << alpha << " " << beta << " " << E << " " << EE << " " << r12 <<endl;
}

void VMC::WriteAccepts(string filename) {
	// Writing the number of accepted steps, assuming Run or Run_NoSaved has been called
	ofstream file;
	file.open("../data/" + filename, ios_base::app);

	file << MCCs << " " << alpha << " " << step_length << " " << accepted <<endl;
}

void VMC::WriteVirial(string filename) {
	ofstream file;
	file.open("../data/" + filename, ios_base::app);	

	double E = ExpectationValues[0]/MCCs;
	double EE = ExpectationValues[1]/MCCs;
	double r12 = ExpectationValues[2]/MCCs;
	double pot = ExpectationValues[3]/MCCs;
	double r12_inverse = ExpectationValues[4]/MCCs;
	file << omega << " " << E << " " << EE << " " << r12 << " " << pot << " " << r12_inverse <<endl;
}