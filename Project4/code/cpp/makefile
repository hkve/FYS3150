model: compile_model run_model
2x2 : 2x2_compile 2x2_run

compile_model:
	g++ -O3 -o IsingModel2D.out IsingModel2D.cpp

run_model: 
	./IsingModel2D.out 

2x2_compile:
	g++ -O3 -o 2x2_lattice.out 2x2_lattice.cpp IsingModel2D.cpp
2x2_run:
	./2x2_lattice.out $s $e $f