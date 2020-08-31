#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "super_general.hpp"
#include "general.hpp"

using namespace std;

void super_general(string method, int max_p);
void general(string method, int max_p);

int main(int argc, char const *argv[])
{
	int max_p;
	string method;

	// Reading method for data and the maximum exponent for matrix dims
	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << " also add methode (sgeneral, general or specific) and max power of n. ex: ./main general 4 \n";
		exit(1);
	}
	else {
		method = argv[1]; // Gets method
		max_p = atoi(argv[2]);  // Gets max power and make int
	}


	if(method == "sgeneral") {
		super_general(method, max_p);
	}
	else if(method == "general") {
		general(method, max_p);
	}
	else if(method == "specific") {
		cout << method <<endl;
	}
	else {
		cout << method << " is not a valid method, enter sgeneral, general or specific" <<endl;
	}

	return 0;
}

void super_general(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		S_general problem(p);
		problem.Initialize();	
		problem.LU();
		problem.Forward_sub();
		problem.Backward_sub();
		
		problem.Write_to_file(method);
		problem.Delete();
	}
}

void general(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		General problem(p); // Call constructor
		problem.Initialize(); // Fill arrays

		auto start = chrono::steady_clock::now();
		
		problem.Forward_sub(); 
		problem.Backward_sub();

		auto end = chrono::steady_clock::now();
		
		cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << " ns" <<endl;

		problem.Write_to_file(method); // Write results
		problem.Delete(); // Free up memory
	}
}