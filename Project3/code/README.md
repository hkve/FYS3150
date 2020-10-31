# Project 3 code
### This folder contains all the programs used to solve project 3
It should contain the following programs:
*main.cpp*
*master.py*
*getInitialConditions.py*
*file_reader.py*
*stability.py*
*escape_velocity.py*
*full_solar_system.py*
*precession.py*

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
- *master.py*  
This program manages *main.cpp* and makes compiling and running the program simpler.

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
  

- *stability.py*  
This program runs tests on our algorithms. It has three main functions:  
  * plot_circular_orbit(dt, T_end, method, N_write)  
  Plots the orbit of the Earth around the Sun for initial conditions which should
  a circular orbit.
    * dt (float): integration stepsize
    * T_end (float): time to simulate forward
    * method (string): either "euler" or "verlet"
    * N_write (int): Number of integration points to write to file
  * plot_error(N_start, N_end, n_tests)  
