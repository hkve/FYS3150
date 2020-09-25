#include "WriteEigs.hpp"

#include <iostream>
#include <cmath>
#include <iomanip>

#include <string>
#include <fstream>

using namespace std;

void writeToFile(string filename, int N, double* eigvals, double** eigvecs, int header_argc, double* header_args) { 
	filename = "data/" + filename + ".dat";

	ofstream outfile (filename, ios_base::app); // Create file

    for (int i=0; i<header_argc; i++) {
        outfile << setw(20) << setprecision(8) << header_args[i];
    }
    outfile << endl;

	for(int i = 0; i < N; i++) {
		outfile << setw(20) << setprecision(8) << eigvals[i];
		for(int j = 0; j < N; j++) {
			outfile << setw(20) << setprecision(8) << eigvecs[i][j];
		}
		outfile << endl;
	}
	outfile.close();
}