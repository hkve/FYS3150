#ifndef WRITEEIGS_HPP
#define WRITEEIGS_HPP

#include <string>

using namespace std;

void writeToFile(string filename, int N, double* eigvals, double** eigvecs, int header_argc, double* header_args);

#endif