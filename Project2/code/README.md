This folder contains all the programs used to solve project 2.
It should contain the following:
JacobiEigSolver.cpp
JacobiEigSolver.hpp
WriteEigs.cpp
WriteEigs.hpp
BucklingBeam.cpp
QuantumOscillator.cpp
file_reader.py
plots.py

- JacobiEigSolver
The central algorithm for this project is the Jacobi Rotation method for
solving eigenvalue problems. This is implemented with the
JacobiEigSolver class which takes an NxN symmetric matrix A as input.
The solver runs until all the diagonal elements are less than a
tolerance, which can be set with the setTolerance method.

- WriteEigs
This program implements a bit of a modular file writing algorithm to
write the eigenvalues and eigenvectors of a solved problem to a
specified filename in the data-folder. Each data entry contains the
paramaters of the specific problem and the number of Jacobi-iterations
necessary to solve it. Since the Quantum problems require a greater
amount of parameters, this function was necessary, as a means to write
the appropriate number of header values in the appropriate places.

BucklingBeam
This implements the specific problem for solving the buckling beam
differential equation. When running the program it requires an argument
which it uses as the N, to set up the NxN matrix representation of the
differential equation. It then uses JacobiEigSolver to solve the
equation and WriteEigs to write these to a data file.

QuantumOscillator
Same as for BucklingBeam, this program implements the specific matrices
for the differential equations relevant to the harmonic oscillator
potential with either one or two electrons. It uses the JacobiEigSolver
class to solve the matrix equaiton and write to file using WriteEigs.
It requires three arguments when run as: number_of_electrons (1 or 2), N
(specifying the size of the matrix to solve), rho_max (the cut off
point for the discretization of rho), and optionally omega_r (required
for the 2-electron system,is set to 1.0 by default)