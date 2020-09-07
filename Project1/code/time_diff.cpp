#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <iomanip> 
#include "LU.hpp"
#include "thomas.hpp"
#include "thomas_singval.hpp"

using namespace std;

double lu(int p);
double thomas(int p);
double thomas_singval(int p);

void write_to_file(double **times, int n_iter, int n_methodes, string filename);

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

		for(int p = 1; p <= max_p_LU; p++) {
			times[i][p-1] = lu(p);
		}

		for(int p = 1; p <= max_p; p++) {
			times[i][max_p_LU+p-1] = thomas(p);
			times[i][max_p_LU+max_p+p-1] = thomas_singval(p);
		}
	}

	write_to_file(times, n_iter, n_methodes, "times");
	return 0;
}


double lu(int p) {
	LU problem(p);
	problem.Initialize();

	auto start = chrono::steady_clock::now();
	problem.Decomp();
	problem.Forward_sub();
	problem.Backward_sub();		
	auto end = chrono::steady_clock::now();
	double time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

	problem.Delete();

	return time;
}

double thomas(int p) {
	Thomas problem(p); // Call constructor
	problem.Initialize(); // Fill arrays

	auto start = chrono::steady_clock::now();
	problem.Forward_sub(); 
	problem.Backward_sub(); 
	auto end = chrono::steady_clock::now();
	double time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();

	problem.Delete(); 
	return time;
}

double thomas_singval(int p) {
	Thomas_singval problem(p);
	problem.Initialize();

	auto start = chrono::steady_clock::now();
	problem.Backward_sub();
	auto end = chrono::steady_clock::now();
	double time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	
	problem.Delete();
	return time; // make milisec
}
	
void write_to_file(double **times, int n_iter, int n_methodes, string filename) {
	filename = "data/" + filename  + ".txt";

	ifstream ifile(filename);
	if(ifile) { // Check if file exists
		remove(filename.c_str()); // In that case, remove it 
	}

	ofstream outfile (filename); // Create file

	for(int i = 0; i < n_iter; i++) {
		for(int j = 0; j < n_methodes; j++) {
			outfile << setprecision(8) << times[i][j] << " ";
		}
		outfile <<endl;
	}
	outfile.close();
}