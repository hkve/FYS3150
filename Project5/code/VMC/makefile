# Only include programs used in the report
all: alpha

# All C++ should be compiled to "compiled" folder
model:
	g++ -o ../compiled/main.exe variationalMC.cpp 

alpha: 
	g++ -o ../compiled/variateAlpha.exe variateAlpha.cpp variationalMC.cpp

sandbox:
	g++ -o ../compiled/sandbox.exe sandbox.cpp variationalMC.cpp