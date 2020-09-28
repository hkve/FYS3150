**This folder contains all the programs used to solve project 2.**  
It should contain the following:  
*JacobiEigSolver.cpp*  
*JacobiEigSolver.hpp*  
*WriteEigs.cpp*  
*WriteEigs.hpp*  
*BucklingBeam.cpp*  
*QuantumOscillator.cpp*  
*file_reader.py*  
*plots.py*  

## C++ programs
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

- BucklingBeam
This implements the specific problem for solving the buckling beam
differential equation. When running the program it requires an argument
which it uses as the N, to set up the NxN matrix representation of the
differential equation. It then uses JacobiEigSolver to solve the
equation and WriteEigs to write these to a data file.

- QuantumOscillator
Same as for BucklingBeam, this program implements the specific matrices
for the differential equations relevant to the harmonic oscillator
potential with either one or two electrons. It uses the JacobiEigSolver
class to solve the matrix equaiton and write to file using WriteEigs.
It requires three arguments when run as: number_of_electrons (1 or 2), N
(specifying the size of the matrix to solve), rho_max (the cut off
point for the discretization of rho), and optionally omega_r (required
for the 2-electron system,is set to 1.0 by default)


## Python programs
- plots.py
This file contains the functions used for all the data representation in
the project. It will check if there are uncompiled C++-files and compile
them, and has functionality to run the files for certain parameters if
the particular data-plot requires it.
The file is run through flags, all which should follow the same '-'
sign. The available flags are:  
    * -h : List the available flags with short descriptions  
    -v : Plot the first couple of eigenvectors of the Buckling Beam
         problem. If the data required is not generated, it will use the
         BucklingBeam program to create a suitable data set.  
    * -c : Plot the number of iterations required for JacobiEigSolver to
         converge on the desired results logarithmically, using the
         BucklingBeam data as reference. A linear fit is also made on the
         data, to determine the convergence rate.  
    * -t : Plot the time performance of our JacobiEigSolver versus the
         Armadillo library's eig_sym-function to find eigenvalues and
         eigenvectors of a symmetric matrix. It will automatically
         generate the appropriate dataset if not already generated.  
    * -r : Plot a 3D plot of the maximum error in the four first
         eigenvalues for the QuantumOscillator for a single electron
         system. The data set is generated automatically if not already
         present.  
    * -q : Plot the specified (arg 2) eigenvector of data from the
         QuantumOscillator of the electron system specified (arg 1).
         If the two-electron system is specified, the remaining arguments
         specify what dataset to generate, before plotting the specifed
         eigenvector. In this case, the data will automatically generate
         using omega_r [0.01, 0.5, 1, 5].  
         args:  
        - no_electrons (str) - either 'one' or 'two', specifying for
                what system to plot.  
        - n (int) - must be greater than 0. Cannot exceed the smallest
                N used for generating the data. Specifies the energy
                level for which to plot the eigenfunction.  
        - N (int) : only for 'two'-electron plotting. Specifies the
                N used for generating the data to plot.  
        - rho_max_list (float) : the last 4 arguments must be a list
                of rho_maxes to use for the omega_rs, in the order as
                above.  
    * -e : Print a LaTeX table of the data from the QuantumOscillator for
         the specified system. It will only inlcude the eigenvalues
         specified in the args, but will not run any program to generate
         the data (yet).  
         args:  
        - no_electrons (str) - either 'one' or 'two', specifying for
                what system to print values.  
        - start_idx (int) - the index value for the first eigenvalue to
                print.  
        - stop_idx (idx) - the index value for the last eigenvalue to
                print.  
