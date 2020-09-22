all: main

main:
	g++ -O3 -o main.exe main.cpp JacobiEigSolver.cpp -larmadillo
	./main.exe  

bb:
	g++ -O3 -o BucklingBeam.exe BucklingBeam.cpp JacobiEigSolver.cpp -larmadillo
	./BucklingBeam.exe