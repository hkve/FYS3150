import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file


def check_init(filename, body_dict, fixedCoM):	# Check if init file exists
	if not filename in os.listdir("initData"): 				
		setInitialConditions(filename, body_dict, fixedCoM)
	else:
		pass

def has_data(filenames):
	if type(filenames) != list:
		filenames = [filenames]

	DATADIR = os.listdir("data")
	for filename in filenames:
		if filename not in DATADIR:
			return False

	return True

def energy(body):
	m = body.m

	R = np.linalg.norm(body.r, axis=0)
	V = np.linalg.norm(body.v, axis=0)
	
	return 0.5*V**2 - 1/R**2

def angular_momentum(body):
	r = body.r
	v = body.v

	L = np.cross(r, v, axis=0)
	L = np.linalg.norm(L, axis=0)

	return L

def plot_circular_orbit(dt=0.0001, T_end=1, method="verlet", N_write=1000):
	"""
	Makes plot of stable Sun/Earth orbit
	Sun at origin no init vel, earth x = 1 AU vx = 2pi * AU/yr

	Args:
		dt: (float) time step in days
		T_end: (int/float) end time in years
		method: (string) "euler" or "verlet"
		N_write: (int) number of points to write to file
	"""
	N = int(T_end/dt)
	
	if N_write == None:
		N_write = N

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	initFilename = "SunEarthStable_init.dat"
	outFilename = "SunEarthStable_" + "_".join([str(method), str(T_end), str(N), str(N_write)])  + ".dat"

	check_init(initFilename, body_dict)

	exists = has_data(outFilename)

	if not exists:
		master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						-out {outFilename} -Nwrite {N_write} -time years {dt} {N}" 
		subprocess.call(master_call.split())

	
	system = read_data_file(outFilename)	

	r = system["Earth"].r

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()

		l = 1.5
		ax.set(xlim=(-l,l), ylim=(-l,l), xlabel="[AU]", ylabel="[AU]")
		ax.axis("equal")
		
		ax.plot(r[0], r[1], c="k", label="Earth") # Earth
		ax.scatter(0,0, c="r", s=150, label="Sun") # Sun
		ax.legend()
		plt.show()
	

def plot_error(dt_start=-1, dt_end=-8, n_tests=30):
	dt = np.logspace(dt_start, dt_end, n_tests, endpoint=True)
	
	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{dt_start}_{dt_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	check_init(initFilename, body_dict, fixedCoM=True)

	exists = has_data(outFilenames)
	i = 0

	if not exists:
		for method in methods:
			for dt_ in dt:
				n_ = int(1/dt_)
				master_call = f"python3 master.py {dt_} {n_} -method {method} -sys initData/{initFilename} -out {outFilenames[i]} -Nwrite {2} -time years --q" 
				subprocess.call(master_call.split())
				i += 1

	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames:
		system = read_data_file(outfile)
		
		rS = system["Sun"].r
		rE = system["Earth"].r

		r0 = np.linalg.norm(rE[:,0]+rS[:,0])
		rN = np.linalg.norm(rE[:,-1]+rS[:,-1])
		error = np.abs((r0-rN)/r0)

		if system["method"] == 0:
			eulerError.append(error)
		if system["method"] == 1:
			verletError.append(error)
	
	with sns.axes_style("darkgrid"):
		fix, ax = plt.subplots()
		ax.set(xscale="log", yscale="log", xlabel="dt", ylabel="Relative error")
		ax.invert_xaxis()
		ax.scatter(dt, eulerError, c="blue", label="Euler")

		ax.scatter(dt, verletError, c="orange", label="Verlet")
		
		# Only do linfit if standard params (as some points have to be excluded)
		if all([val is plot_error.__defaults__[i] for i, val in enumerate([dt_start,dt_end,n_tests])]): 
			sE, eE = 3, len(dt) # Euler points to linfit (start, end)
			sV, eV = 0, 11      # Verlet points to linfit (start, end)
			
			# Linfit for Euler
			slopeE, constE, r_valueE, p_valueE, std_errE = \
			stats.linregress(np.log10(dt[sE:eE]), np.log10(eulerError[sE:eE]))

			# Linfit for Verlet
			slopeV, constV, r_valueV, p_valueV, std_errV = \
			stats.linregress(np.log10(dt[sV:eV]), np.log10(verletError[sV:eV]))

			# Plot Euler and Verlet linfit
			ax.plot(dt[sE:eE], 10**constE * dt[sE:eE]**slopeE, c="red",\
					label=f"Slope Euler = {slopeE:.3f}$\pm${std_errE:.3f}", markersize=3)
			ax.plot(dt[sV:eV], 10**constV * dt[sV:eV]**slopeV, c="k" ,\
					label=f"Slope Verlet = {slopeV:.3f}$\pm${std_errV:.3f}", markersize=3)
		

		ax.legend()
		plt.show()
	

def plot_energy(N=int(1e7), T_end = 50, N_write=10000):
	dt = T_end/N
	
	methods = ["euler", "verlet"]
	initFilename = "SunEarthStable_init.dat"

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}

	check_init(initFilename, body_dict, fixedCoM=True)

	outFilenames = []
	for method in methods:
		outFilenames.append(f"SunEarthStable_{method}_{T_end}_{N}_{N_write}.dat")

	exists = has_data(outFilenames)

	i = 0
	if not exists:
		for method in methods:
			master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
							 -out {outFilenames[i]} -Nwrite {N_write} -time years {dt} {N}" 
			subprocess.call(master_call.split())
			i += 1


	with sns.axes_style("darkgrid"):
		t = np.linspace(0, T_end, N_write)
		fig, ax = plt.subplots()
		ax.set(yscale="log", ylim=(1e-4, 1))
		c = ["k","r"]
		for i, method in enumerate(methods):
			system = read_data_file(outFilenames[i])
			E = energy(system["Earth"])
			L = angular_momentum(system["Earth"])
			
			E = np.abs((E-E[0])/E[0])
			#L = np.abs((L-L[0])/L[0])
			ax.plot(t, E, label=r"$E_{tot} $ " + method.capitalize(), c=c[i])
			ax.plot(t, L, label=r"$L$ " + method.capitalize())
	ax.legend()
	plt.show()

def plot_time(N_start=2, N_end=8, n_tests=30):
	N = np.logspace(N_start, N_end, n_tests, endpoint=True, dtype=int)
	
	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{N_start}_{N_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	check_init(initFilename, body_dict, fixedCoM=True)

	exists = has_data(outFilenames)
	i = 0

	if not exists:
		for method in methods:
			for n_ in N:
				dt_ = 1/n_
				master_call = f"python3 master.py {dt_} {n_} -method {method} -sys initData/{initFilename} -out {outFilenames[i]} -Nwrite {2} -time years --q" 
				subprocess.call(master_call.split())
				i += 1

	eulerTime = []
	verletTime = []
	for outfile in outFilenames:
		system = read_data_file(outfile)
		
		if system["method"] == 0:
			eulerTime.append(system["time"])
		if system["method"] == 1:
			verletTime.append(system["time"])


	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(xscale="log", yscale="log")
		ax.set_xlabel(xlabel="N", fontsize=13)
		ax.set_ylabel(ylabel="Time [s]", fontsize=13)
		ax.scatter(N, eulerTime, label="Euler")
		ax.scatter(N, verletTime, label="Verlet")
	
		# Only linfit if standard params (as some points have to be excluded)
		if all([val is plot_time.__defaults__[i] for i, val in enumerate([N_start,N_end,n_tests])]):
			sE, eE = 15, len(N)
			sV, eV = 15, len(N)

			slopeE, constE, r_valueE, p_valueE, std_errE = \
			stats.linregress(np.log10(N[sE:eE]), np.log10(eulerTime[sE:eE]))

			# Linfit for Verlet
			slopeV, constV, r_valueV, p_valueV, std_errV = \
			stats.linregress(np.log10(N[sV:eV]), np.log10(verletTime[sV:eV]))

			# Plot Euler and Verlet linfit
			ax.plot(N[sE:eE], 10**constE * N[sE:eE]**slopeE, c="red",\
					label=f"Slope Euler = {slopeE:.3f}$\pm${std_errE:.3f}", markersize=3)
			ax.plot(N[sV:eV], 10**constV * N[sV:eV]**slopeV, c="k" ,\
					label=f"Slope Verlet = {slopeV:.3f}$\pm${std_errV:.3f}", markersize=3)
		
	ax.legend(fontsize=13)
	plt.show()

plot_time()