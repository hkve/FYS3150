import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import getInitialCondition, setInitialConditions
from file_reader import read_data_file

def plot_Sun_Earth_Jupiter(dt=0.001, T = 12, jupiter_scale = 1):
	rj = 5.2 
	vj = 2*np.pi/np.sqrt(rj)
	body_dict = {
		"Sun": [0,0,0,0,0,0],
		"Earth": [1,0,0,0,2*np.pi,0],
		"Jupiter": [-rj,0,0,0,vj,0]
	}
	scaled_mass = {"Jupiter": jupiter_scale}

	initFilename = f"SunEarthJupiter_stable_{jupiter_scale}m.dat"
	outFilename = f"SunEarthJupiter_stable_{T}_{jupiter_scale}m.dat"
	setInitialConditions(initFilename, body_dict, fixedCoM=False, scaled_mass=scaled_mass)

	N = int(T/dt)
	N_write = 1000

	if True:
		master_call = f"python3 master.py -method verlet -sys initData/{initFilename} \
							-out {outFilename} -Nwrite {N_write} -time years --fixSun {dt} {N}" 
		subprocess.call(master_call.split())

	system = read_data_file(outFilename)
	
	rE, rJ = system["Earth"].r, system["Jupiter"].r
	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.scatter(0,0, c="r", label="Sun")
		ax.plot(rE[0], rE[1], label="Earth")
		ax.plot(rJ[0], rJ[1], label="Jupiter")
		ax.set_xlabel("x [AU]", fontsize=13)
		ax.set_ylabel("y [AU]", fontsize=13)
		ax.axis("equal")
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=3, fancybox=True, shadow=True, fontsize=13)
	plt.show()
	
def plot_radial_distance(dt=0.001, T=12, jupiter_scale=1000):
	rj = 5.2 
	vj = 2*np.pi/np.sqrt(rj)
	body_dict = {
		"Sun": [0,0,0,0,0,0],
		"Earth": [1,0,0,0,2*np.pi,0]
	}

	scale_j_mass = [1, 10, 100, 1000]

	N = int(T/dt)
	N_write = 10000
	t = np.linspace(0, T, N_write)
	initFilenameNoJ = "No_jupiter.dat"
	outFilenameNoJ = "No_jupiter_data.dat"
	
	# First run without jupiter
	setInitialConditions(initFilenameNoJ, body_dict, fixedCoM=False)
	master_call = f"python3 master.py -method verlet -sys initData/{initFilenameNoJ} \
					-out {outFilenameNoJ} -Nwrite {N_write} -time years --fixSun {dt} {N}" 
	subprocess.call(master_call.split())

	system = read_data_file(outFilenameNoJ)
	rE_no_j = np.linalg.norm(system["Earth"].r, axis=0)
	rE = []
	# Add Jupiter to the mix
	body_dict["Jupiter"] = [-rj,0,0,0,-vj,0]
	identifier = "Jupiter_inf" 
	for scale in scale_j_mass:
		filename = identifier + str(scale) + ".dat"
		scaled_mass = {"Jupiter": scale}
		setInitialConditions(filename, body_dict, scaled_mass=scaled_mass)

		master_call = f"python3 master.py -method verlet -sys initData/{filename} \
					 -out {filename} -Nwrite {N_write} -time years --fixSun --q {dt} {N}" 
		subprocess.call(master_call.split())

		system = read_data_file(filename)
		R = np.linalg.norm(system["Earth"].r, axis=0)
		rE.append(R)

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		for R, lab in zip(rE, scale_j_mass):
			if lab == 1:
				ax.plot(t, rE_no_j-R, label="$M_J$")
			else:
				ax.plot(t, rE_no_j-R, label=f"{lab}$M_j$")

		ax.set_xlabel("Time [years]", fontsize=13)
		ax.set_ylabel("Deviation from circular orbit [AU]", fontsize=13)
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
          ncol=4, fancybox=True, shadow=True)

	plt.show()


plot_Sun_Earth_Jupiter(jupiter_scale=1000)