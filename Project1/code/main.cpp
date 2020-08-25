#include "problem.hpp"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char const *argv[])
{

	int max_n;
	string filename;

	if (argc <= 1) {
		cout << "bad usage: " << argv[0] << ", also add filename and a max power of n. ex: ./main outfile 4" << "\n";
		exit(1);
	}
	else {
		filename = argv[1];
		max_n = atoi(argv[2]); 
	}

	for(int i = 1; i <= max_n; i++) {
		Problem Problem1(i);	

		Problem1.Initialize();
		Problem1.LU();
		Problem1.Forward_substitution();
		Problem1.Backward_substitution();
		Problem1.Write_to_file(filename);

		//Problem1.Print_problem();
		//Problem1.Print_LU();
		Problem1.Print_sol();
		Problem1.Delete();	
	}

	return 0;
}