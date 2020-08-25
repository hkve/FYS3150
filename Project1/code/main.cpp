#include "problem.hpp"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char const *argv[])
{

	int max_n;
	string filename;

	// Reading filename for data and the maximum exponent for matrix dims
	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << ", also add filename and a max power of n. ex: ./main outfile 4" << "\n";
		exit(1);
	}
	else {
		filename = argv[1];
		max_n = atoi(argv[2]);  // Make int
	}

	for(int i = 1; i <= max_n; i++) {
		Problem Problem1(i); // Make object	

		Problem1.Initialize(); // Filling array and seting up vectors
		Problem1.LU(); // Preform LU decomposition
		Problem1.Forward_substitution(); // Solving Ly = b_tilde
		Problem1.Backward_substitution(); // Solving LUv = b_tilde => Uv = y
		Problem1.Write_to_file(filename); // Write results to file

		// Some print functions
		//Problem1.Print_problem();
		//Problem1.Print_LU();
		//Problem1.Print_sol();
		Problem1.Delete();	// Free up memory (i hope) for the next exponent
	}

	return 0;
}