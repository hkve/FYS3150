#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{   
    // tror kommentarene under er feil, evt fjern de senere
    
	/*args:
    T (double): Temperature in kT/J
    MCCs_start (double): Cycle start 
    MCCs_end (double): Cycle end
    dMCCs (int): Number of MCCs between each MCCs calculated between start and end
    initSpin (int): 0 or 1. The orientation of the init spins. 0 for random, 1 for all positively aligned
    */
    
	double T = atof(argv[1]);
    int MCCs = atof(argv[2]);
    int initSpin = atoi(argv[3]);
	string filename = (string) argv[4];//"../data/20x20Values.out";
	//string filename = "../data/20x20Values.out";
    IsingModel* problem = new IsingModel(20,MCCs,1,T);
    problem->Initialize(initSpin);
    problem->Solve();
    problem->writeFinalExpValues(filename);
    // int MCCs_start = atof(argv[2]);
    // int MCCs_end = atof(argv[3]);
    // int dMCCs = atoi(argv[4]);
    //string outfile = (string)argv[6];


	

	return 0;
}