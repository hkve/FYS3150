import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file


def check_init(filename, body_dict):	# Check if init file exists
	if not filename in os.listdir("initData"): 				
		setInitialConditions(filename, body_dict)
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

def plot_circular_orbit(dt=0.001, T_end=1, method="verlet", N_write=1000):
	"""
	Makes plot of stable Sun/Earth orbit
	Sun at origin no init vel, earth x = 1 AU vx = 2pi * AU/yr

	Args:
		dt: (float) time step in days
		T_end: (int/float) end time in years
		method: (string) "euler" or "verlet"
		N_write: (int) number of points to write to file
	"""
	d2y = 365  	# Day to year conversion
	N = int(T_end*d2y/dt)
	
	if N_write == None:
		N_write = N

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi/d2y,0]}
	
	initFilename = "SunEarthStable_init.dat"
	outFilename = "SunEarthStable_" + "_".join([str(method), str(T_end), str(N), str(N_write)])  + ".dat"

	check_init(initFilename, body_dict)

	exists = has_data(outFilename)

	if exists == False:
		master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						-out {outFilename} -dpts {N_write} {dt} {N}" 
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
	

def plot_error(N_start=3, N_end=7, n_tests=30):
	d2y = 365.25

	log_start, log_end = np.log10(N_start), np.log10(N_end)
	N = np.logspace(log_start, log_end, n_tests)
	N = (10**N).astype(int)
	
	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{N_start}_{N_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi/d2y,0]}
	

	check_init(initFilename, body_dict)
	exists = has_data(outFilenames)
	i = 0
	if exists == False:
		for method in methods:
			for n in N:
				dt = d2y/n

				master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						    	-out {outFilenames[i]} -dpts {1} {dt} {n}" 
				subprocess.call(master_call.split())
				i += 1

	
	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames:
		system = read_data_file(outfile)
		rx = system["Earth"].r[0]
		error = np.abs(rx[0]-rx[-1])

		if system["method"] == 0:
			eulerError.append(error)
		if system["method"] == 1:
			verletError.append(error)
	
	with sns.axes_style("darkgrid"):
		fix, ax = plt.subplots()
		ax.set(xscale="log", yscale="log", xlabel="N", ylabel="Relative error")
		ax.scatter(N, eulerError, label="Euler")
		ax.scatter(N, verletError, label="Verlet")
		ax.legend()
		plt.show()
	

def plot_energy(N=int(1e7), T_end = 50, N_write=10000):
	d2y = 365
	dt = T_end*d2y/N
	
	methods = ["euler", "verlet"]
	initFilename = "SunEarthStable_init.dat"

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi/d2y,0]}

	check_init(initFilename, body_dict)

	outFilenames = []
	for method in methods:
		outFilenames.append(f"SunEarthStable_{method}_{T_end}_{N}_{N_write}.dat")

	exists = has_data(outFilenames)

	i = 0
	if exists == False:
		for method in methods:
			master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
							 -out {outFilenames[i]} -dpts {N_write} {dt} {N}" 
			subprocess.call(master_call.split())
			i += 1


	with sns.axes_style("darkgrid"):
		t = np.linspace(0, T_end, N_write+1)
		fig, ax = plt.subplots()
		ax.set(yscale="log", ylim=(1e-4, 1))
		c = ["k","r"]
		for i, method in enumerate(methods):
			system = read_data_file(outFilenames[i])
			E = energy(system["Earth"])
			E = np.abs((E-E[0])/E[0])
			ax.plot(t[:-1], E[:-1], label=r"$E_{tot} $ " + method.capitalize(), c=c[i])
	ax.legend()
	plt.show()

plot_error()