main: compile_main run_main
bb: compile_bb run_bb

# For testing, migth remove
compile_main:
	g++ -O3 -o main.exe main.cpp JacobiEigSolver.cpp
	  
run_main:
	./main.exe $N

# Buckling Beam 
compile_bb:
	g++ -O3 -o BucklingBeam.exe BucklingBeam.cpp JacobiEigSolver.cpp
run_bb:
	./BucklingBeam.exe $N
