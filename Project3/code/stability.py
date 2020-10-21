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

	if exists == False:
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
	

def plot_error(N_start=1, N_end=9, n_tests=10):
	log_start, log_end = np.log10(N_start), np.log10(N_end)
	N = np.logspace(log_start, log_end, n_tests, endpoint=True)
	N = (10**N).astype(int)
	

	methods = ["euler", "verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		for i in range(n_tests):
			outFilenames.append(f"SunEarthStable_{method}_{N_start}_{N_end}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	setInitialConditions(initFilename, body_dict)

	#check_init(initFilename, body_dict)
	exists = has_data(outFilenames)
	i = 0

	if not True:#exists == False:
		for method in methods:
			for n in N:
				dt = 1/n
				master_call = f"python3 master.py {dt} {n} -method {method} -sys initData/{initFilename} -out {outFilenames[i]} -Nwrite {2} -time years --q" 
				subprocess.call(master_call.split())
				i += 1

	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames:
		try:
			system = read_data_file(outfile)
			
			r = system["Earth"].r
			#print(r, "r", r[0,-1], r[1,-1])
			error = np.abs(r[0,0] - r[0,-1])#np.abs(1-np.sqrt((r[0,-1])**2 + r[1,-1]**2))
			print(error)
			

			if system["method"] == 0:
				# print("euler error: ", rx[0], rx[-1])
				eulerError.append(error)
			if system["method"] == 1:
				# print("verlet error: ", rx[0], rx[-1])
				verletError.append(error)
		except:
			verletError.append(1)
	
	with sns.axes_style("darkgrid"):
		fix, ax = plt.subplots()
		ax.set(xscale="log", yscale="log", xlabel="N", ylabel="Relative error")
		
		ax.scatter(N, eulerError, label="Euler")
		ax.scatter(N, verletError, label="Verlet")
		ax.legend()
		plt.show()
	

def plot_energy(N=int(1e7), T_end = 50, N_write=10000):
	dt = T_end/N
	
	methods = ["euler", "verlet"]
	initFilename = "SunEarthStable_init.dat"

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}

	check_init(initFilename, body_dict)

	outFilenames = []
	for method in methods:
		outFilenames.append(f"SunEarthStable_{method}_{T_end}_{N}_{N_write}.dat")

	exists = has_data(outFilenames)

	i = 0
	if exists == False:
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

plot_error()