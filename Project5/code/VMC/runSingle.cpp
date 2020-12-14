#include "variationalMC.hpp"
#include "trialFunctions.cpp"

int main(int argc, char *argv[]) {
	if(argc < 6) {
		cerr << "Bad usage enter:" <<endl;
		cerr << "Trial wavefunction: T1 or T2" <<endl;
		cerr << "Local energy function: E0, E1 or E2" <<endl;
		cerr << "log10(MCCs)" <<endl;
		cerr << "omega" <<endl;
		cerr << "Outfilename" <<endl;
		exit(1);
	}

	string trialFunction = argv[1];
	string Local_energy = argv[2];
	int MCCs = pow(10, atof(argv[3]));
	double omega = atof(argv[4]);
	string outfilename = argv[5];

	func psi;
	func EL;

	double alpha;
	double beta;

	// Choose trial function
	if(trialFunction == "T1") {
		psi = psi_T1;
	}
	else if(trialFunction == "T2") {
		psi = psi_T2;
	}
	else {
		cerr << trialFunction << " is not a recognised trial function. Use T1 or T2" <<endl;
		exit(1);
	}

	// Choose local energy
	if(Local_energy == "E0") {
		EL = EL_1;
		alpha = 1; beta = 0;
	}
	else if(Local_energy == "E1") {
		EL = EL_1_Coulomb;
		alpha = 0.9; beta = 0;
	}
	else if(Local_energy == "E2") {
		EL = EL_2;
		alpha = 0.994212; beta = 0.286315;
		//alpha = 0.95; beta = 0.34;
	}
	else {
		cerr << Local_energy << " is not a recognised local energy. Use E0, E1 or E2" <<endl;	
		exit(1);
	}

	VMC* problem = new VMC(psi, EL, MCCs, 0, omega, alpha, beta);
	problem->Run_NoSave();
	problem->WriteVariational(outfilename);
}