# Only include programs used in the report
all: alpha

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