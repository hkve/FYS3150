#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{   
    // tror kommentarene under er feil, evt fjern de senere
    
	/*args:
    T (double): Temperature in kT/J
    MCCs (int): Number of cycles to run 
    initSpin (int): 0 or 1. The orientation of the init spins. 0 for random, 1 for all positively aligned
    filename (string) Name of file to save in data folder
    */
    
	double T = atof(argv[1]);
    int MCCs = atof(argv[2]);
    int initSpin = atoi(argv[3]);
	string filename = (string) argv[4];
    IsingModel* problem = new IsingModel(20,MCCs,1,T);
    problem->Initialize(initSpin);
    problem->Solve_ExpValuesMeanwhile(filename); // This saves ExpValues at every cycle


	

	return 0;
}