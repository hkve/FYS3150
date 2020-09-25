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

# Quantum Oscillator
compile_qo:
	g++ -O3 -o QuantumOscillator.exe QuantumOscillator.cpp QuantumOscillator.cpp WriteEigs.cpp
# Run for 1 electron
run_qo1:
	./Quantum.exe 1 $N $r
# Run for 2 electrons
run_qo1:
	./Quantum.exe 2 $N $r $o