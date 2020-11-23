# Project 4 code
### This folder contains all the programs used to solve project 4
The cpp directory should contain the following files:
*IsingModel2D.cpp*,
*IsingModel2D.hpp*,
*IsingModel2DParalell.cpp*,
*IsingModel2DParalell.hpp*,
*2x2_lattice.cpp*,
*20x20_lattice.cpp*,
*20x20_energy_count.cpp*,
*acceptedFlips.cpp*,
*parallel.cpp*,
*timeComparison.cpp*,
*plot_spins.cpp*

The python directory should contain the following files:

*main.py*,
*plotExpValues2x2.py*,
*plotMCS2x2.py*,
*plot20x20.py*,
*plotAcceptedFlips.py*,
*plotEnergyHistogram.py*,
*plotPhaseTransition.py*,
*plotCVphasetrans.py*,

## Usage
### Installs
The following packages must be installed to run **all** of the code. Some will run with missing installs.

- `numpy`
- `matplotlib`
- `argparse`
- `colour`
- `scipy`


### main.py
Located in the **python** directory, *main.py* serves as the controll panel for this code base. This file takes arguments from the command line and can be used to run the ising model for different lattice sizes both in series and in parallel as well as recreate most of the data and plots shown in the report. **NB:** make sure to add --compile as a flag if you want to run a specific program/plot and have not compiled. If you want to compile everything simply navigate to the **python** directory and write.

```console
foo@bar:~$ python3 main.py compile
```

To see what kind of simulations you can run simply type and what kind of plots you can create simply type.

```console
foo@bar:~$ python3 main.py -h
```
Some plots shown in the report are not available through **main.py** due to taking a very long time to compute. If you want to create these plots navigate to the **python** directory and open the desired .py program. Preform a function call of main with sim = True. **NB:** this will take a lot of time.

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
