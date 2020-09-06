#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include "LU.hpp"
#include "thomas.hpp"
#include "special.hpp"

using namespace std;

void lu(string method, int max_p);
void thomas(string method, int max_p);
void special(string method, int max_p);

int main(int argc, char const *argv[])
{
	int max_p;
	string method;

	// Reading method for data and the maximum exponent for matrix dims
	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << " also add method (LU, thomas or special) and max power of n. ex: ./main thomas 4 \n";
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
	else if(method == "special") {
		special(method, max_p);
	}
	else {
		cout << method << " is not a valid method, enter LU, Thomas or Thomas_singular" <<endl;
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

		auto start = chrono::steady_clock::now();
		
		problem.Forward_sub(); 
		problem.Backward_sub();

		auto end = chrono::steady_clock::now();
		
		cout << "p = " << p << " t = ";
		cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << " ns" <<endl;

		problem.Write_to_file(method); // Write results
		problem.Delete(); // Free up memory
	}
}

void special(string method, int max_p) {
	for(int p = 1; p <= max_p; p++) {
		Special problem(p);
		problem.Initialize();

		auto start = chrono::steady_clock::now();
		
		problem.Backward_sub();

		auto end = chrono::steady_clock::now();
		cout << "p = " << p << " t = ";
		cout << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << " ns" <<endl;

		problem.Write_to_file(method);
		problem.Delete();
	}
}