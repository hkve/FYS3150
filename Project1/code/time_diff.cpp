#include <iostream>
#include <fstream>
#include <chrono>
#include "LU.hpp"
#include "thomas.hpp"
#include "thomas_singval.hpp"

using namespace std;

double lu(int p);
//double thomas(int p);
//double thomas_singval(int p);

int main(int argc, char const *argv[])
{
	int n_iter;
	int n_methodes;
	int max_p;
	int max_p_LU;

	// Reading method for data and the maximum exponent for matrix dims
	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << " also add number for iterations and max power of n (LU will stop at 3). ex: ./main 500 7 \n";
		exit(1);
	}
	else {
		n_iter = atoi(argv[1]);
		max_p = atoi(argv[2]);  // Gets max power and make int
	}

	if(max_p > 3) {
		n_methodes = 3 + 2*max_p;
		max_p_LU = 3;
	}
	else {
		n_methodes = 3*max_p;
		max_p_LU = max_p;
	}
	
	double** times = new double*[n_iter];
	for(int i = 0; i < n_iter; i++) {
		times[i] = new double[n_methodes];

		for(int p = 1; p < max_p_LU; p++) {
			times[i][p] = lu(p);
			cout << times[i][p] <<endl;
		}
	}
	return 0;
}


double lu(int p) {
	LU problem(p);
	problem.Initialize();

	auto current_time = std::chrono::system_clock::now();
	problem.Decomp();
	problem.Forward_sub();
	problem.Backward_sub();		
	auto duration_in_seconds = std::chrono::duration<double>(current_time.time_since_epoch());
	problem.Delete();

	double num_seconds = duration_in_seconds.count();
	return num_seconds;
}
/*
double thomas(int p) {
	Thomas problem(p); // Call constructor
	problem.Initialize(); // Fill arrays
	problem.Forward_sub(); 
	problem.Backward_sub(); 
	problem.Delete(); 
}

double thomas_singval(int p) {
	Thomas_singval problem(p);
	problem.Initialize();
	problem.Backward_sub();
	problem.Delete();
}
*/