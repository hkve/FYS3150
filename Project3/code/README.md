# Project 3 code
### This folder contains all the programs used to solve project 3
It should contain the following programs:
*main.cpp*,
*master.py*,
*plotter.py*,
*getInitialConditions.py*,
*file_reader.py*,
*stability.py*,
*escape_velocity.py*,
*jupiter_influence.py*,
*full_system.py*,
*modified_gravity.py*,
*precession.py*

## Usage
### Installs
The following packages must be installed to run master.py and plotter.py. 

- `numpy`
- `pandas`
- `matplotlib`
- `colour`
- `seaborn`
- `astroquery`


### master.py
To run simulations using *main.cpp* we recommend to use the *master.py*. This file takes arguments from the command line and runs *main.cpp*
with those parameters. **NB:** make sure to add --compile as a flag if main.exe is missing. The parameters are:
- required
    * **dt** (int/float): log10 of step length for the simulation in AU/yr
    * **N** (int/float): log10 of the number of integration steps to use
- optional
    * **h, --help** Shows help message (equal to this information) and exits.
    * **-method** Integration method to use. Must be "euler" or "verlet". Default: "verlet".
    * **-beta** Change the inverse proportionality of gravity. Must be in range [2,3]. Default = 2 (normal, newtonian gravity)
    * **-sys** Name of init file in initdData/ dir. Default: sys.dat
    * **-out** Name of file in data/ dir where simulation results are stored. Default: sys.out
    * **-Nwrite** Numbers of data points to be stored/written. Defaults to N. NB, not logarithmic!
    * **--GR** Do simulation with general relativity correction term.
    * **--fixSun**  Do simulation with sun fixed at 0,0,0. (Sun must be first element in init file!)
    * **--compile** Compile main.cpp to "main.exe" before running
    * **--q** Run quietly

### Example
Provided that the file *my_sys.dat* and *my_data.out* exsits in initData and data respectively, this call compiles *main.cpp* to *main.exe*, runs *main.exe* with step length h=1e-4 and N=1e5, writing 1000 points with the sun fixed at the center 
```console
foo@bar:~$ python3 master.py -4 5 -sys my_sys.dat -out my_data.out -Nwrite 1000 --fixSun --compile
```
### plotter.py
To recreate the figures shown in the report, we recomend to use *plotter.py*. This file also takes arguments from the commandline, corresponding to different plots. All programs run the required simulation and reads the data dynamically. If you want to change parameters, go into the separate plot files (look in *plotter.py* to find in what files the different plotting code is located) and change default parameters. The parameters set equals the ones used in the report. The plots available are:
- optional
  * **-h, --help**         show this help message and exit
  * **-circularOrbit**     Plot stable circulat orbit of earth, with sun kept fixed
  * **-ellipticalOrbits**  Plot elliptical orbits of earth for different inital velocities, with sun kept fixed
  * **-benchmark**         Preforms a benchmarking of Euler and velocity Verlet and plots the result
  * **-error**             Check energy conservation for Euler and velocity Verlet and plots the result
  * **-sunEarthJupiter**   Plot Sun-Earth-Jupiter system with the mass of jupiter scaled by a facotr of 1000
  * **-radialDistance**    Plot the radial deviation from Earths orbit without Jupiter for different masses of jupiter
  * **-escapeVelocity**    Plot a visual representation of escape velocity for different inital velocities
  * **-fullSystem**        Plot our own real system
  * **-beta10yr**          Path of planet for beta=2.92 over 10 years
  * **-betaPosEnergy**     Plot path and energy of several beta-sytems

This call will plot the orbit of Earth around the sun with initial conditions creating a circular orbit over 1 year, and with a modfied inverse proportionality of beta = 2.92 (figure 1 and 13 in the report respectively)  
### Example
```console
foo@bar:~$ python3 plotter.py -circularOrbit -beta10yr
```

# Brief explanation of central programs
## C++ programs
- *main.cpp*
This program contains the struct Body and the class System.
  * Body:  
    Contains all the information relevant to a body in the
    N-body system. It contains three dimensional position and
    velocity as well as the mass and the UUID.
  * System:  
    This class contains all the methods and algorithms necessary to
    solve the N-body system. It reads from an init-file and creates
    a Body for every line in the init-file, with the appropriate
    initial conditions specified in the init-file. The solve method
    solves the N-body system with the stepsize specified, making
    N-1 steps forward with the method specified. The results are
    written to a specified filename.  
args:  
    * initfile (string): the file location of the init-file which
    specifies the system and the initial conditions.  
    * outfile (string): the filename where the data written to file
    will be stored.  
    * Nwrite (int): The number of datapoints to be written to file.
    The endpoints are always included, and the points between are
    space as evenly as possible.  
    * dt (double): Argument given as an exponent of 10. Serves as the
    stepsize used to solve.  
    * N (int): Argument given as an exponent of 10. the system will be
    solved for N-1 steps forward.  
    * method (int): Specifies the numerical solving method. 0 for the
    forward Euler method, 1 for velocity Verlet.  
    * beta (double): Specifies the proportionality of the gravitational
    force. Default is normal gravity, beta = 2.0.  
    * GR (bool): Specifies whether the general relativity correction term
    to Newtonian gravity should be added to the force. Default is false.  
    * timeFormat (string): Either "years" or "days". Specifies whether the
    initial conditions measures time in years or days.  
    * q (bool): If true, the program will not output anything to terminal.

## Python programs
- *getInitialConditions.py*  
This program links with the *astroquery*-API to the JPL Horizons database. It can
get the desired data from the database, and translate it into the init-files used
by our solver. It can also set up such init-files using dictionaries with the
desired initial conditions. It has two main functions:
  * getInitialCondition(filename, bodies, date, fixedCoM)  
  Gets the initial conditions of the specified bodies on the specified date and writes
  an init-file with the specified filename.
    * filename (string): filename for the init-file to be written. It will automatically
    be placed in the initData folder.
    * bodies (list(string)): A list of the names of the bodies to be included 
    (i.e. "Sun", "Earth") in the init file. It will default to include all bodies
    * date (string): YYYY-MM-DD format. The date for the initial conditions data to
    be pulled. Defaults to today.
    * fixedCoM (bool): Whether to adjust the initial conditions such that the center
    of mass will remain fixed at the origin. Defaults to false.
    * scaled_mass (dict): A dictionary where the keys are the ordinary names of the bodies and value is 
    a number to scale the mass by. Defualt to None 
  * setInitialConditions(filename, body_dict, fixedCoM)  
  Writes the init-file with the specified filename to the initData-folder with the bodies
  and initial conditions specfied in body_dict.
    * filename (string): Filename given to the init-file.
    * body_dict (dict): A dictionary where the keys are the ordinary names of the bodies
    and the values are lists like [x, y, z, vx, vy, vz].
    * fixedCoM (bool): Whether to adjust the initial conditions such that the center
    of mass will remain fixed at the origin. Defaults to false.
    * scaled_mass (dict): A dictionary where the keys are the ordinary names of the bodies and value is 
    a number to scale the mass by. Defualt to None 
    
- *file_reader.py*  
This programs simply reads the datafiles created by main.cpp. The data is returned as a dictionary containing all 
the relevant information about that datafile. 
   * Body (class): Each body read from the file is stored as an instance of this class. It contains properties such as the name of the body,
   its mass, UUID, as well as the position  and velocity vectors.
   * read_data_file(filename):
   Reads a data file and stores each body as an instance of the body class, containg all the information conserning that body. 
   It also stores information about that specific simulation such as the method used, step length, number of integration steps,
   time used to preform the simulation and how many points that are written to the actual data file
     * filename (string): The name of the file to read from, stored in the data directory 
