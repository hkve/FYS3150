import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file


def check_init(filename, body_dict, fixedCoM=False):	# Check if init file exists
	if not filename in os.listdir("initData"): 				
		setInitialConditions(filename, body_dict, fixedCoM=fixedCoM)
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
	
	return 0.5*V**2 , -1/R

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

	check_init(initFilename, body_dict, fixedCoM=False)
	setInitialConditions(initFilename, body_dict, fixedCoM=False)

	exists = has_data(outFilename)

	if not exists:
		master_call = f"python3 master.py -method {method} -sys initData/{initFilename} \
						-out {outFilename} -Nwrite {N_write} -time years {dt} {N}" 
		subprocess.call(master_call.split())

	
	system = read_data_file(outFilename)	
	r = system["Earth"].r

	Ek, Ep = energy(system["Earth"])
	Ek_std, Ep_std = np.std(Ek), np.std(Ep)
	Ek_mean, Ep_mean = np.mean(Ek), np.mean(Ep)
	
	if system["method"] == 0:
		method = "euler"
	if system["method"] == 1:
		method = "verlet"

	print(f"Method = {method}, dt = {dt:E}, N={N:E}, T_end = {T_end}")
	print(f"Ek initial = {Ek[0]:.4f}, Ep initial = {Ep[0]:.4f}")
	print(f"Ek (mean) = {Ek_mean:.4f}, std = {Ek_std:.4E}, relative = {(Ek_std/Ek_mean):.4E}")
	print(f"Ek (mean) = {Ep_mean:.4f}, std = {Ep_std:.4E}, relative = {(Ep_std/Ep_mean):.4E}")

	T = np.linspace(0, T_end, N_write, endpoint=True)
	plt.plot(T, Ek)
	plt.show()
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



def plot_elliptical_orbits(dt=0.0001, T_end=2000, n_v = 4, method="verlet", N_write=100000):
	N = int(T_end/dt)
	if N_write == None:
		N_write = N

	VX = np.linspace(np.pi, 5*np.pi/2, n_v)

	filenames = [f"SunEarthEllip_{n_v}_{i}_{T_end}.dat" for i in range(n_v)]
	
	for i in range(n_v):
		body_dict = {"Sun": [0,0,0,0,0,0],
				    "Earth": [1,0,0,0,VX[i],0]}
		check_init(filenames[i], body_dict, fixedCoM=False)
		
	exists = has_data(filenames)

	if not exists:
		for i in range(n_v):
			master_call = f"python3 master.py -method {method} -sys initData/{filenames[i]} \
							-out {filenames[i]} -Nwrite {N_write} -time years --q {dt} {N}" 
			subprocess.call(master_call.split())

	labs = ["$\pi$", "$3\pi/2$", "$2\pi$", "$5\pi/2$"]

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set_xlabel("x [AU]", fontsize=13)
		ax.set_ylabel("y [AU]", fontsize=13)
		ax.axis("equal")
		for i in range(n_v):
			system = read_data_file(filenames[i])
			rE = system["Earth"].r
			ax.plot(rE[0], rE[1], label=f"$v_x=${labs[i]} AU/yr")

		ax.scatter(0,0, c="r")
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.18),
          ncol=2, fancybox=True, shadow=True, fontsize=13)
		plt.show()


	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		T = np.linspace(0, T_end, N_write, endpoint=True)
		for i in range(n_v):
			system = read_data_file(filenames[i])
			L = angular_momentum(system["Earth"])
			ax.plot(T, L)
			
	plt.show()


def plot_energy(N=int(1e4), T_end = 50, N_write=10000):
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
		ax.set(yscale="log", ylim=(1e-14, 1))
		c = ["k","r"]
		for i, method in enumerate(methods):
			system = read_data_file(outFilenames[i])
			Ep, Ek = energy(system["Earth"])
			E = Ek+Ep
			L = angular_momentum(system["Earth"])
			
			E = np.abs((E-E[0])/E[0])
			L = np.abs((L-L[0])/L[0])
			print(L)
			ax.plot(t, E, label=r"$E_{tot} $ " + method.capitalize(), c=c[i])
			ax.plot(t, L, label=r"$L$ " + method.capitalize())
	ax.legend()
	plt.show()



def plot_time(N_start=2, N_end=8, n_tests=50):
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
	tot = 2*n_tests
	if not exists:
		for method in methods:
			for n_ in N:
				dt_ = 1/n_
				master_call = f"python3 master.py {dt_} {n_} -method {method} -sys initData/{initFilename} -out {outFilenames[i]} -Nwrite {2} -time years --q" 
				subprocess.call(master_call.split())
				i += 1
				print(f"Done {method}, dt = {n_}, {(i*100/tot):.2f}%")
	eulerTime = []
	verletTime = []
	for outfile in outFilenames:
		system = read_data_file(outfile)
		
		if system["method"] == 0:
			eulerTime.append(system["time"])
		if system["method"] == 1:
			verletTime.append(system["time"])

	eulerTime = np.array(eulerTime)
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
		if all([val is plot_time.__defaults__[i] for i, val in enumerate([N_start,N_end,n_tests])]):
			sE, eE = 25, len(N)
			sV, eV = 25, len(N)

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


def plot_Etot_error(dt_start=3, dt_end=7, n_tests=20):
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
	if not exists:
		for method in methods:
			for N_, dt_ in zip(N,dt):
				master_call = f"python3 master.py {dt_} {N_} -method {method} -sys initData/{initFilename} \
							   -out {outFilenames[i]} -Nwrite {2} -time years --log --fixSun --q" 
				subprocess.call(master_call.split())
				i += 1
				print(f"{method}, log10(dt)={dt_:.2f}, log10(N)={N_:.2f}, {(i*100/tot):.1f}%")

	N = 10**N
	dt = 10**dt
	systems = []
	eulerError = []
	verletError = []

	for outfile in outFilenames:
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

