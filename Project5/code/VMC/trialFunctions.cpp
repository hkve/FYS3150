/*
Here you can define all the trial wavefuntions you want
Use the absolute square of the WF
REMEMBER: in the file you want to use them add #include "trialFunctions.cpp"
NB: Don't add trialFunctions.cpp to compile list, cpp does this auto
*/

double psi_T1(double *r, double omega, double alpha, double beta) {
	// Add the trial wavefunction
	double r1_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double r2_squared = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
	
	return exp(-alpha*omega*(r1_squared+r2_squared));
}

double psi_T2(double *r, double omega, double alpha, double beta) {
	double r12 = (r[0]-r[3])*(r[0]-r[3]) +
				 (r[1]-r[4])*(r[1]-r[4]) +
				 (r[2]-r[5])*(r[2]-r[5]);
	r12 = sqrt(r12);

	return psi_T1(r, omega, alpha, beta)*exp(r12/(1+beta*r12));	
}

double EL_1(double *r, double omega, double alpha, double beta) {
	// Add the local energy, ie (1/|PSI_T>)H|PSI_T>
	double r1_squared = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double r2_squared = r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
	
	return 0.5*omega*omega*(r1_squared+r2_squared)*(1-alpha*alpha) + 3*alpha * omega;
}

double EL_1_Columb(double *r, double omega, double alpha, double beta) {
	double r12 = (r[0]-r[3])*(r[0]-r[3]) +
				 (r[1]-r[4])*(r[1]-r[4]) +
				 (r[2]-r[5])*(r[2]-r[5]);
	return EL_1(r, omega, alpha, beta) + 1/sqrt(r12);
}

double EL_2(double *r, double omega, double alpha, double beta) {
	double r12 = (r[0]-r[3])*(r[0]-r[3]) +
				 (r[1]-r[4])*(r[1]-r[4]) +
				 (r[2]-r[5])*(r[2]-r[5]);
	r12 = sqrt(r12);
	double betaTerm = 1/(1+beta*r12);
	double betaTermSquared = 0.5*betaTerm*betaTerm;

	double extraAddition = alpha*omega*r12 - betaTermSquared - 2/r12 + 2*betaTerm;
	// Think calling EL_1 here is faster, no need to calculate r12 twice
	return EL_1(r, omega, alpha, beta) + 1/r12 + betaTermSquared*extraAddition;	
}