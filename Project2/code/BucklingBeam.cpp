#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
	if(argc <= 1) {
		cout << "Bad usage, enter matrix dim NxN, example ./BucklingBeam 20";
		exit();
	}
	return 0;
}