#include "IsingModel2D.hpp"

int main(int argc, char const *argv[])
{
	if(argc <= 2) {
		cerr << "bad usage, enter Tstart, Tend and output filename. (Example ./2times2.out 1 4 expValues.out";
		exit(1);
	}
	cout << argv[0] << " " << argv[1] << " " << argv[2];
	return 0;
}