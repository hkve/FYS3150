/*
Here you can define all the trial wavefuntions you want
REMEMBER: in the file you want to use them add #include "trialFunctions.cpp"
NB: Don't add trialFunctions.cpp to compile list, cpp does this auto
*/

double psi_T1(double *r, double omega, double alpha) {
	// Add the trial wavefunction
	cout << "i am the trial wavefunction" <<endl;
	return 0;
}

double EL_1(double *r, double omega, double alpha) {
	// Add the local energy, ie (1/|PSI_T>)H|PSI_T>
	cout << "and i am the local energy!" <<endl;
	return 0;
}