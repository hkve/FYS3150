compile_model:
	g++ -O3 -o IsingModel2D.out IsingModel2D.cpp

2x2:
	g++ -O3 -o 2x2_lattice.out 2x2_lattice.cpp IsingModel2D.cpp

acceptedFlips:
	g++ -O3 -o acceptedFlips.out acceptedFlips.cpp IsingModel2D.cpp

20x20:
	g++ -O3 -o 20x20_lattice.out 20x20_lattice.cpp IsingModel2D.cpp

energyCount:
	g++ -O3 -o 20x20_energy_count.out 20x20_energy_count.cpp IsingModel2D.cpp

parallel:
	g++ -O3 -o parallel.out parallel.cpp IsingModel2DParalell.cpp -fopenmp