compile_model:
	g++ -o compiled/main.exe variationalMC.cpp 

compile_alpha: 
	g++ -o compiled/variateAlphaNoCoulomb.exe variateAlphaNoCoulomb.cpp variationalMC.cpp

compile_sandbox:
	g++ -o compiled/sandbox.exe sandbox.cpp variationalMC.cpp