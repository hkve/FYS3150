# Only include programs used in the report
all: alpha optimalStep beta betaGrid virial single

# All C++ should be compiled to "compiled" folder
model:
	g++ -O3 -o ../compiled/main.exe variationalMC.cpp 

alpha: 
	g++ -O3 -o ../compiled/variateAlpha.exe variateAlpha.cpp variationalMC.cpp

sandbox:
	g++ -O3 -o ../compiled/sandbox.exe sandbox.cpp variationalMC.cpp

optimalStep:
	g++ -O3 -o ../compiled/optimalStep.exe optimalStep.cpp variationalMC.cpp

beta:
	g++ -O3 -o ../compiled/variateBeta.exe variateBeta.cpp variationalMC.cpp

betaGrid:
	g++ -O3 -o ../compiled/variateBeta_grid.exe variateBeta_grid.cpp variationalMC.cpp 

parallel:
	g++ -O3 -o ../compiled/variateBeta_grid_parallel.exe variateBeta_grid_parallel.cpp variationalMC.cpp -fopenmp

virial:
	g++ -O3 -o ../compiled/virialTheorem.exe virialTheorem.cpp variationalMC.cpp

single:
	g++ -O3 -o ../compiled/runSingle.exe runSingle.cpp variationalMC.cpp