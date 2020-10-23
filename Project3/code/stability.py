import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file


def check_init(filename, body_dict):	# Check if init file exists
	if not filename in os.listdir("initData"): 				
		setInitialConditions(filename, body_dict, fixedCoM=True)
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
	N = int(1e5)
	dt = 1e-5
	N_write = 2
	
	if N_write == None:
		N_write = N

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	initFilename = "SunEarthStable_init.dat"
	#initFilename = "FSS_init.dat"
	outFilename = "SunEarthStable_" + "_".join([str(method), str(T_end), str(N), str(N_write)])  + ".dat"

	#check_init(initFilename, body_dict)
	setInitialConditions(initFilename, body_dict, fixedCoM=True)

	#exists = has_data(outFilename)

	if True:#not exists:
		master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						-out {outFilename} -Nwrite {N_write} -time years {dt} {N}" 
		subprocess.call(master_call.split())

	
	system = read_data_file(outFilename)	
	print(system)
	r = system["Earth"].r
	rs = system["Sun"].r
	#rm = system["Mercury"].r
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()

		l = 1.5
		ax.set(xlim=(-l,l), ylim=(-l,l), xlabel="[AU]", ylabel="[AU]")
		ax.axis("equal")
		
		ax.plot(r[0], r[1], "-o", c="r", label="Earth") # Earth
		ax.plot(rs[0], rs[1], "-o", c="k", label="Sun") # Earth

		#ax.plot(rm[0], rm[1], c="b", label="Mercury") # Earth

		ax.legend()
		plt.show()
	

def plot_error(N_start=1, N_end=3, n_tests=20):
	log_start, log_end = np.log10(N_start), np.log10(N_end)
	N = np.linspace(N_start, N_end, n_tests, endpoint=True) 
	#N = np.array([1,2,3,4,5,6,7,8], dtype=np.float64)
	#N_start = N[0]
	#N_end = N[-1]
	dt = -N
	#N = (10**N).astype(int)
	

	methods = ["verlet"]

	initFilename = "SunEarthStable_init.dat"

	outFilenames = []
	for method in methods:
		# for i in range(n_tests):
		for i in range(len(N)):
			outFilenames.append(f"SunEarthStable_{method}_{int(N_start)}_{int(N_end)}_{i+1}.dat") 

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,0,2*np.pi,0]}
	
	setInitialConditions(initFilename, body_dict, fixedCoM=True)

	#check_init(initFilename, body_dict)
	#exists = has_data(outFilenames)
	i = 0

	if not False:#not exists:
		for method in methods:
			for n, dt_ in zip(N,dt):
				#dt = 1/n
				master_call = f"python3 master.py {dt_} {n} -method {method} -sys initData/{initFilename} -out {outFilenames[i]} -Nwrite {2} -time years --log" 
				subprocess.call(master_call.split())
				i += 1
	N = 10**N
	dt = 10**dt
	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames:
	
		system = read_data_file(outfile)
		def f(x,y):
			return np.linalg.norm(x-y, axis=0)
		r = system["Earth"].r
		rs = system["Sun"].r
		r0 = r[:,0]
		r1 = r[:,-1]
		rs0 = rs[:,0]
		rs1 = rs[:,-1]

		#print(r, "r", r[0,-1], r[1,-1])
		pos0 = np.sqrt((r[0,0]-rs[0,0])**2+ (r[1,0]-rs[1,0])**2)    #np.abs(1-np.sqrt((r[0,-1])**2 + r[1,-1]**2))
		pos1 = np.sqrt((r[0,-1]-rs[0,-1])**2+ (r[1,-1]-rs[1,-1])**2)
		#pos0 = np.abs(r[0,0]-rs[0,0])
		#pos1 = np.abs(r[0,-1]-rs[0,-1])
		
		error = np.abs((pos1-pos0)/pos0)
		error = np.abs(f(rs1, r1)-f(rs0, r0))
		#print(error)
		#print("")
		
		# 17 18
		if system["method"] == 0:
			# print("euler error: ", rx[0], rx[-1])
			eulerError.append(error)
		if system["method"] == 1:
			# print("verlet error: ", rx[0], rx[-1])
			verletError.append(error)

	
	with sns.axes_style("darkgrid"):
		fix, ax = plt.subplots()
		ax.set(xscale="log", yscale="log", xlabel="N", ylabel="Relative error")
		verlim = -7
		vs = 0
		#fitEuler = np.polyfit(np.log10(N[5:]), np.log10(eulerError[5:]), deg = 1)
		fitVerlet = np.polyfit(np.log10(N[vs:verlim]), np.log10(verletError[vs:verlim]), deg = 1)
		#print("euler fit: ", fitEuler)
		print("verlet fit: ", fitVerlet)
		#ax.scatter(N, eulerError, label="Euler")
		ax.scatter(N, verletError, label="Verlet")
		ax.plot(N[vs:verlim], N[vs:verlim]**(fitVerlet[0]) * 10**fitVerlet[1])
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

# plot_error()
plot_circular_orbit()