import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file
from master import simulate

def check_init(filename, body_dict):
	"""
	Check if init file exsists, in case not make it
		Args:
			Filename: (string) What the init file should be named
			Body_dict: (dictionary) Holding the initial conditions for each body
	"""
	if not filename in os.listdir("initData"): 				
		setInitialConditions(filename, body_dict)
	else:
		pass

def has_data(filenames): 
	"""
	Check if data files exsists
		Args: 
			filenames (string/list of strings) Name of the files that should be there
		Returns:
			True: if all files are present
			False: if one or more files are missing
	"""
	if type(filenames) != list:
		filenames = [filenames]

	DATADIR = os.listdir("data")
	for filename in filenames:
		if filename not in DATADIR:
			return False

	return True

def energy(body): 
	"""
	Calculates the kinetic and potential energy of a body
		Args: 
			Body: (Body) what body to calculate on
		Returns:
			Ek: Kinetick energy of body in [ME(AU/yr)^2]
			EP: Potential energy of body in [ME(AU/yr)^2]
	"""
	m = body.m

	R = np.linalg.norm(body.r, axis=0)
	V = np.linalg.norm(body.v, axis=0)
	 
	return 0.5*V**2 , -4*np.pi**2/R # [ME(AU/yr)^2]

def angular_momentum(body):
	"""
	Calculates the angular momentum of a body
		Args:
			Body: (Body) what bodt to calculate on
		Returns:
			L: (float array) Array of the length of angular momentum vector over time
	"""
	r = body.r
	v = body.v

	L = np.cross(r, v, axis=0)
	L = np.linalg.norm(L, axis=0)

	return L


def circularOrbit(dt=0.0001, T_end=1, method="verlet", N_write=1000):
	"""
	Makes plot of stable Sun/Earth orbit
	Sun at origin no init vel, earth x = 1 AU vx = 2pi * AU/yr

	Args:
		dt: (float) time step in years
		T_end: (int/float) end time in years
		method: (string) "euler" or "verlet"
		N_write: (int) number of points to write to file
	"""
	N = int(T_end/dt)
	if N_write == None:
		N_write = N

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]} # Initial conditions
	
	initFilename = "SunEarthStable_init.dat"
	outFilename = "SunEarthStable_" + "_".join([str(method), str(T_end), str(N), str(N_write)])  + ".dat" # Make filenames

	check_init(initFilename, body_dict) 
	setInitialConditions(initFilename, body_dict)

	exists = has_data(outFilename)
	N = np.log10(N)
	dt = np.log10(dt)
	if not exists: # If the datafiles are missing, make them
		simulate(N=N, dt = dt, method=method, Nwrite=N_write, sys=initFilename, out=outFilename, fixSun=True, quiet=True)

	
	system = read_data_file(outFilename) # Reads the data
	r = system["Earth"].r

	Ek, Ep = energy(system["Earth"]) # Calculate energy
	Ek_std, Ep_std = np.std(Ek), np.std(Ep)
	Ek_mean, Ep_mean = np.mean(Ek), np.mean(Ep)
	
	if system["method"] == 0:
		method = "euler"
	if system["method"] == 1:
		method = "verlet"

	print(f"Method = {method}, dt = {dt:E}, N={N:E}, T_end = {T_end}")
	print(f"Ek initial = {Ek[0]:.4f}, Ep initial = {Ep[0]:.4f}")
	print(f"Ek (mean) = {Ek_mean:.4f}, std = {Ek_std:.4E}, relative = {(Ek_std/Ek_mean):.4E}")
	print(f"Ep (mean) = {Ep_mean:.4f}, std = {Ep_std:.4E}, relative = {(-Ep_std/Ep_mean):.4E}")

	T = np.linspace(0, T_end, N_write, endpoint=True)
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		l = 1.5
		ax.set(xlim=(-l,l), ylim=(-l,l))
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.axis("equal")
		ax.set_xlabel("x [AU]", fontsize=13)
		ax.set_ylabel("y [AU]", fontsize=13)

		ax.scatter(0,0, c="r", label="Sun")
		ax.plot(r[0], r[1], c="b", label="Earth")

		ax.legend(fontsize=13)
		plt.show()

def ellipticalOrbits(dt=0.0001, T_end=10, n_v = 4, method="verlet", N_write=10000):
	"""
	Makes plot of different elliptical orbits
		Args: 
			dt: (float) Step length in yr
			T_end: (int/float) Total time to run the simulation
			n_v: (int) Number of orbits in the range v_y = [pi, 5pi/2]
			method: (string) Verlet or Euler
			N_write: (int) How many points to write
	"""
	N = int(T_end/dt)
	if N_write == None:
		N_write = N

	VY = np.linspace(np.pi, 5*np.pi/2, n_v, endpoint=True) # Array holding different initial velocities in the y-direction

	filenames = [f"SunEarthEllip_{n_v}_{i}_{T_end}.dat" for i in range(n_v)] # Storing filenames
	
	for i in range(n_v):
		body_dict = {"Sun": [0,0,0,0,0,0],
				    "Earth": [1,0,0,0,VY[i],0]}
		check_init(filenames[i], body_dict)
		
	exists = has_data(filenames)
	dt = np.log10(dt)
	N = np.log10(N)
	if not exists:
		for i in range(n_v):
			simulate(N=N, dt = dt, method=method, Nwrite=N_write, sys=filenames[i], out=filenames[i], fixSun=True, quiet=True)

	labs = ["$\pi$", "$3\pi/2$", "$2\pi$", "$5\pi/2$"] # Only true if standard parameters are used

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set_xlabel("x [AU]", fontsize=14)
		ax.set_ylabel("y [AU]", fontsize=14)
		ax.axis("equal")
		for i in range(n_v):
			system = read_data_file(filenames[i])
			rE = system["Earth"].r
			ax.plot(rE[0], rE[1], label=f"$v_y=${labs[i]} AU/yr")

		ax.scatter(0,0, c="r")
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.22),
          ncol=2, fancybox=True, shadow=True, fontsize=13)
		plt.show()


	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		T = np.linspace(0, T_end, N_write, endpoint=True)
		ax.set_xlabel("Time [yr]", fontsize=14)
		ax.set_ylabel("Relative error of |$\ell$|", fontsize=14)
		for i in range(n_v):
			system = read_data_file(filenames[i])
			L = angular_momentum(system["Earth"])
			L = np.abs((L[0]-L)/L[0])
			ax.plot(T, L, label=f"$v_y=${labs[i]}")
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18),
          ncol=2, fancybox=True, shadow=True, fontsize=13)
	plt.show()


def benchmark(N_start=2, N_end=7, n_tests=50):
	"""
	Preform a timing of the two algorithms and make a log-log plot
		Args:
			N_start: (int) log10 of the number of integration steps to start at
			N_end: (int) log10 of the number of integration steps to end at
			n_tests: (int) number of test between N_start and N_end
	"""
	N = np.log10(np.logspace(N_start, N_end, n_tests, endpoint=True, dtype=int))

	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{N_start}_{N_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	check_init(initFilename, body_dict)

	exists = has_data(outFilenames)

	i = 0
	tot = 2*n_tests 
	if not exists: # If datafiles are missing, make them
		for method in methods:
			for n_ in N:
				dt_ = -n_
				simulate(N=n_, dt = dt_, method=method, Nwrite=2, sys=initFilename, out=outFilenames[i], fixSun=True, quiet=True)
				i += 1
				print(f"Done {method}, dt = {n_}, {(i*100/tot):.2f}%") # It takes a while... good to know where your at
	eulerTime = []
	verletTime = []

	for outfile in outFilenames: # Sort through files and place the time in euler and verlet
		system = read_data_file(outfile)
		
		if system["method"] == 0:
			eulerTime.append(system["time"])
		if system["method"] == 1:
			verletTime.append(system["time"])

	eulerTime = np.array(eulerTime) # For easier handling
	verletTime = np.array(verletTime)
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(xscale="log", yscale="log")
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.set_xlabel(xlabel="N", fontsize=15)
		ax.set_ylabel(ylabel="Time taken [s]", fontsize=15)
		ax.scatter(N, eulerTime, label="Euler")
		ax.scatter(N, verletTime, label="Verlet")

		# Only linfit if standard params (as some points have to be excluded)
		if all([val is benchmark.__defaults__[i] for i, val in enumerate([N_start,N_end,n_tests])]):
			sE, eE = 25, len(N) # What points to fit Euler
			sV, eV = 25, len(N) # What point to fit Verlet

			slopeE, constE, r_valueE, p_valueE, std_errE = stats.linregress(np.log10(N[sE:eE]), np.log10(eulerTime[sE:eE]))

			# Linfit for Verlet
			slopeV, constV, r_valueV, p_valueV, std_errV = \
			stats.linregress(np.log10(N[sV:eV]), np.log10(verletTime[sV:eV]))

			# Plot Euler and Verlet linfit
			ax.plot(N[sE:eE], 10**constE * N[sE:eE]**slopeE, c="red",\
					label=f"Slope Euler = {slopeE:.3f}$\pm${std_errE:.3f} \n $R^2$ = {(r_valueE**2):.4f}", markersize=3)
			ax.plot(N[sV:eV], 10**constV * N[sV:eV]**slopeV, c="k" ,\
					label=f"Slope Verlet = {slopeV:.3f}$\pm${std_errV:.3f} \n $R^2$ = {(r_valueV**2):.4f}", markersize=3)
			
	ax.legend(fontsize=13)
	plt.show()


def error(dt_start=3, dt_end=7, n_tests=20):
	"""
	Making a plot of the relative error in the total energy of the two algorithms
		Args:
			dt_start: (int) -log10 of what dt to start at
			dt_end: (int) -log10 of what dt to end at 
			n_tests: (int) number of tests to preform between dt_start and dt_end
	"""
	dt = np.linspace(dt_start, dt_end, n_tests, endpoint=True)*-1 
	N = -dt	

	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{dt_start}_{dt_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	setInitialConditions(initFilename, body_dict)
	check_init(initFilename, body_dict)
	exists = has_data(outFilenames)
	i = 0
	tot = 2*n_tests
	if not exists: # If files are missing, create them
		for method in methods:
			for N_, dt_ in zip(N,dt):
				simulate(N=N_, dt = dt_, method=method, Nwrite=2, sys=initFilename, out=outFilenames[i], fixSun=True, quiet=True)
				i += 1
				print(f"{method}, log10(dt)={dt_:.2f}, log10(N)={N_:.2f}, {(i*100/tot):.1f}%") # Again, migth take a while

	N = 10**N
	dt = 10**dt
	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames: # Sort through files and place the error in euler and verlet
		system = read_data_file(outfile)
		Ek, Ep = energy(system["Earth"])
		
		E = Ek+Ep
		error = np.abs((E[0]-E[-1])/E[0])
	
		if system["method"] == 0:
			eulerError.append(error)
		if system["method"] == 1:
			verletError.append(error)
	
	with sns.axes_style("darkgrid"):
		fix, ax = plt.subplots()
		ax.invert_xaxis()
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.set_xlabel("$h$ [yr]", fontsize=14)
		ax.set_ylabel("Relative error of $E_{tot}$", fontsize=14)
		ax.set(xscale="log", yscale="log")
		ax.scatter(dt, eulerError, label="$E_{tot}$ Euler")
		ax.scatter(dt, verletError, label="$E_{tot}$ Verlet")
		slope, const, r_value, p_value, std_err = stats.linregress(np.log10(dt), np.log10(eulerError))
		ax.plot(dt, 10**const * dt**slope, c="k" ,linestyle="dashed",\
				label=f"s = {slope:.3f}$\pm${std_err:.3f} \n$R^2$ = {(r_value**2):.5f}")
		
	ax.legend(fontsize=13, loc=6)
	plt.show()
