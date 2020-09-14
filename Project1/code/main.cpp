/*
This file preforms the calculations for alle the 3 different methodes. 
From the commandline it takes an algorithm (LU, Thomas, or Thomas_singval)
and a maximum power max_p. The chosen algorithm will preform the calculations
for 10^1 to 10^max_p. 

The results will be written to one file for each power. The filename
follows Name_of_algorithm + power + .txt
*/

#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "LU.hpp"
#include "thomas.hpp"
#include "thomas_singval.hpp"

using namespace std;

void lu(string method, int max_p);
void thomas(string method, int max_p);
void thomas_singval(string method, int max_p);

int main(int argc, char const *argv[])
{
	int max_p;
	string method;

	// Reading method for data and the maximum exponent for matrix dims
	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << " also add method (LU, Thomas or Thomas_singval) and max power of n. ex: ./main thomas 4 \n";
		exit(1);
	}
	else {
		method = argv[1]; // Gets method
		max_p = atoi(argv[2]);  // Gets max power and make int
	}


	if(method == "LU") {
		lu(method, max_p);
	}
	else if(method == "Thomas") {
		thomas(method, max_p);
	}
	else if(method == "Thomas_singval") {
		thomas_singval(method, max_p);
	}
	else {
		cout << method << " is not a valid method, enter LU, Thomas or Thomas_singval" <<endl;
	}

	return 0;
}

void lu(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		LU problem(p);
		problem.Initialize();	
		problem.Decomp();
		problem.Forward_sub();
		problem.Backward_sub();
		
		problem.Write_to_file(method);
		problem.Delete();
	}
}

void thomas(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		Thomas problem(p); // Call constructor
		problem.Initialize(); // Fill arrays

		problem.Forward_sub(); 
		problem.Backward_sub();

		problem.Write_to_file(method); // Write results
		problem.Delete(); // Free up memory
	}
}

void thomas_singval(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		Thomas_singval problem(p);
		problem.Initialize();		
		
		problem.Backward_sub();

		problem.Write_to_file(method);
		problem.Delete();
	}
}