# Project 5 code
### This folder contains all the programs used to solve project 5
The VMC directory should contain the following files:
- *variationalMC.hpp*
- *variationalMC.cpp*
- *trialFunctions.cpp*
- *runSingle.cpp*
- *optimalStep.cpp*  
- *variateAlpha.cpp*
- *variateBeta.cpp*                
- *variateBeta_grid.cpp*
- *variateBeta_grid_parallel.cpp*
- *virialTheorem.cpp*

The plotting directory should contain the following files:

- *main.py*               
- *plotStability.py* 
- *plotEvarEvsAlpha.py*
- *optimalStep.py*
- *optimalStepFitter.py*
- *plotAcceptancerate.py* 
- *plotVirial.py*
- *plot3DVarBeta.py*      
- *plotVarBetaParallel.py* 
- *tables.py*


## Usage
### Installs
The following packages must be installed to run **all** of the code. Some will run with missing installs.

- `numpy`
- `matplotlib`
- `seaborn`
- `colour`
- `scipy`


### main.py
Located in the **plotting** directory, *main.py* serves as the controll panel for this code base. This file takes arguments from the command line and can be used to run the VMC for different wave function and local energies as well as recreate most of the data and plots shown in the report. **NB:** Make sure to compile before you run a plot/simulation. Simply navigate to the **plotting** directory and write.

```console
foo@bar:~$ python3 main.py compile
```

To see what kind of simulations you can run simply type and what kind of plots you can create simply type.

```console
foo@bar:~$ python3 main.py -h
```

### Examples
If you want to compile and simulate a 2x2 lattice from T0=1 to T1=8 and output the data to data/foo.dat run
```console
foo@bar:~$ python3 main.py 2x2 -T1 8 -filename foo.dat --compile
```
If you want to simulate a 40x40 lattice from T0 = 1 to T1 = 4 with dT = 0.5 and output the data to data/foo.dat, provided that you have already compiled  
```console
foo@bar:~$ python3 main.py parallel -T0 1 -T1 4 -dT 0.5 -filename foo.dat
```
If you want to create a plot of the expectation values of a 2x2 lattice for different temperatures as shown in the report, run
```console
foo@bar:~$ python3 main.py  plot-exp-2x2 --sim --compile
```
