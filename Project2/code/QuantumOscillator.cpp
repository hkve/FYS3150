#include <iostream>
#include <cmath>
#include <string>
#include "JacobiEigSolver.hpp"

#define PI 3.14159265359

using namespace std;

void setA(double** A, int N, double rho_max, int no_electrons, double omega_r);

int main(int argc, char** argv)
{
	bool time;
	if(argc <= 3) {
		cout << "Bad usage: enter number of electrons [1,2], number of integration points [1, inf), integration endpoint [1, inf) and optionally oscillator frequency [0,inf). Example ./QuantumOscillator two 20 5 12";
		exit(1);
	}

    int NUMBER_OF_ELECTRONS;
    if (atoi(argv[1]) == 1) {
        NUMBER_OF_ELECTRONS = 1;
    } else if (atoi(argv[1]) == 2) {
        NUMBER_OF_ELECTRONS = 2;
    } else {
        cout << "Bad usage: the first argument must be the number of electrons in the systen, either 1 or 2.";
        exit(1);
    }
	
	int	N = atoi(argv[2]);

    double rho_max = atof(argv[3]);

    double omega_r;
    if (argc > 3) {
        omega_r = atof(argv[4]);
    } else {
        omega_r = 1;
    }
    if (NUMBER_OF_ELECTRONS==2) {
        cout << "omega_r was not given, set to 1.0 by default.";
    }


	double** A = new double*[N];
	for(int i = 0; i < N; i++) {
		A[i] = new double[N];
	}
	setA(A, N, rho_max, NUMBER_OF_ELECTRONS, omega_r);

	JacobiEigSolver* problem = new JacobiEigSolver(A, N);

	problem -> Solve();

	string filename = NUMBER_OF_ELECTRONS==1 ?  "QuantumOscillator_one" : "QuantumOscillator_two";
	problem -> writeToFile(filename);
	
	return 0;
}

void setA(double** A, int N, double rho_max, int no_electrons = 1, double omega_r = 1.0) {
    double h, hh, a, d, rho_i;
	h = rho_max / (N-1.0);
    hh = h*h;

	a = -1 / hh;
	d = 2 / hh;
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < N; j++) {
			if(i==j) {
                rho_i = i*i * hh;
				A[i][j] = d + omega_r*rho_i;
                if (no_electrons == 2){
                    A[i][j] += 1/rho_i;
                }
			} else if (abs(i-j) == 1) {
				A[i][j] = a;
			} else {
				A[i][j] = 0.0;
			}
		}
	}
}