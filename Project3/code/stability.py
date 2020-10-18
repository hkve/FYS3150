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
	
	K = 0.5*V**2
	P = 1/R

	return K,P

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

	print(r[0])
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()

		l = 1.5
		ax.set(xlim=(-l,l), ylim=(-l,l), xlabel="[AU]", ylabel="[AU]")
		ax.axis("equal")
		
		ax.plot(r[0], r[1], c="k", label="Earth") # Earth
		ax.scatter(0,0, c="r", s=150, label="Sun") # Sun
		ax.legend()
		plt.show()
	

def plot_error(dt_start=10, dt_end=0.0001, n_tests=10):
	d2y = 365

	log_start, log_end = np.log10(dt_start), np.log10(dt_end)
	DT = np.logspace(log_start, log_end, n_tests)

	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{log_start:.0f}_{log_end:.0f}_{i+1}") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi/d2y,0]}
	
	check_init(initFilename, body_dict)
	exists = has_data(outFilenames)
	i = 0
	if exists == False:
		for method in methods:
			for dt in DT:
				N = int(365/dt)
				
				master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						    	-out {outFilenames[i]} -dpts {1} {dt} {N}" 
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
		
	print(DT)
	fix, ax = plt.subplots()
	ax.set(xscale="log", yscale="log")
	ax.scatter(DT, eulerError, label="euler")
	ax.scatter(DT, verletError, label="verlet")
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
		outFilenames.append(f"SunEarthStable_{method}_{T_end}_{np.log10(N):.0f}_{N_write}.dat")

	exists = has_data(outFilenames)

	i = 0
	if exists == False:
		for method in methods:
			master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
							 -out {outFilenames[i]} -dpts {N_write} {dt} {N}" 
			subprocess.call(master_call.split())
			i += 1

	fig, ax = plt.subplots()
	ax.set(yscale="log")
	c = ["k","r"]

	for i, method in enumerate(methods):
		system = read_data_file(outFilenames[i])
		K, P = energy(system["Earth"])
		rel_K = np.abs((K-K[0])/K[0])
		rel_P = np.abs((P-P[0])/P[0])
		ax.plot(rel_K, label=f"K {method}", c=c[i])
		ax.plot(rel_P, label=f"P {method}", c=c[i], linestyle="dashed")
	ax.legend()
	plt.show()

plot_energy(T_end=50, N=int(1e7))