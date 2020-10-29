import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import getInitialCondition
from file_reader import read_data_file

def plot_Sun_Earth_Jupiter(dt=0.01, T = 15, jupiter_scale = 1):
	T *= 365
	N = int(T/dt)
	N_write = 10000

	initFilename = f"SEJ_{jupiter_scale}.dat"
	outFilename = f"SEJ_{T}_{jupiter_scale}_{N_write}.dat"
	bodies = ["Sun", "Earth", "Jupiter"]

	scaled_mass = {"Jupiter": jupiter_scale}
	if not initFilename in os.listdir("initData"):
		getInitialCondition(initFilename, bodies, fixedCoM=False, scaled_mass=scaled_mass)

	if not outFilename in os.listdir("data"):
		master_call = f"python3 master.py -method verlet -sys initData/{initFilename} \
							-out {outFilename} -Nwrite {N_write} -time days --fixSun {dt} {N}" 
		subprocess.call(master_call.split())

	system = read_data_file(outFilename)
	rS = system["Sun"].r
	rE = system["Earth"].r-rS
	rJ = system["Jupiter"].r-rS

	l = np.max(rJ) + np.max(rJ)/10

	bodynames = ["Earth", "Jupiter"]
	bodycolors = ["#0EEB58", "orange"]

	fig = plt.figure()  
	fig.set_facecolor('black')
	ax = fig.add_subplot(111, projection='3d')  
	ax.set_facecolor('black') 
	ax.grid(False) 
	ax.w_xaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
	ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
	ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.0))
	ax.set(xlim=(-l,l), ylim=(-l,l), zlim=(-l,l))
	ax.scatter(0,0, c="yellow", label="Sun")
	ax.plot(rE[0,1:], rE[1,1:], rE[2,1:], c="#0EEB58", label="Earth")
	ax.plot(rJ[0], rJ[1], rJ[2], c="orange", label="Jupiter")
	ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
             ncol=3, fancybox=True, shadow=True, fontsize=15)
	plt.show()

def plot_radial_distance(dt=0.01, T=12):
	scale_j_mass = [1, 10, 100, 1000]

	T *= 365
	N = int(T/dt)
	N_write = 10000
	initFilenameNoJ = "No_jupiter.dat"
	outFilenameNoJ = f"No_jupiter_{T}_{-np.log10(dt)}.dat"
	
	if not initFilenameNoJ in os.listdir("initData"):
		getInitialCondition(initFilenameNoJ, ["Sun", "Earth"])
	if not outFilenameNoJ in os.listdir("data"):
		master_call = f"python3 master.py -method verlet -sys initData/{initFilenameNoJ} \
					-out {outFilenameNoJ} -Nwrite {N_write} -time days --fixSun --q {dt} {N}" 
		subprocess.call(master_call.split())

	system = read_data_file(outFilenameNoJ)
	rE_no_j = np.linalg.norm(system["Earth"].r, axis=0)
	
	initFilenames = [f"SEJ_{scale}.dat" for scale in scale_j_mass]
	outFilenames = [f"SEJ_{T}_{scale}_{N_write}.dat" for scale in scale_j_mass]

	rE = []
	for i in range(len(scale_j_mass)): 
		if not initFilenames[i] in os.listdir("initData"):
			getInitialCondition(initFilenames[i], ["Sun", "Earth", "Jupiter"], scaled_mass={"Jupiter": scale_j_mass[i]})

		if not outFilenames[i] in os.listdir("data"):
			master_call = f"python3 master.py -method verlet -sys initData/{initFilenames[i]} \
						-out {outFilenames[i]} -Nwrite {N_write} -time days --fixSun --q {dt} {N}" 
			subprocess.call(master_call.split())

		system = read_data_file(outFilenames[i])
		R = np.linalg.norm(system["Earth"].r, axis=0)
		rE.append(R)

	t = np.linspace(0, T/365, N_write, endpoint=True)
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		for R, lab in zip(rE, scale_j_mass):
			diff = rE_no_j-R
			if lab == 1:
				ax.plot(t, diff, label="$M_J$")
			else:
				ax.plot(t, diff, label=f"{lab}$M_J$")

			diff = np.abs(diff)
			avg_diff = np.mean(diff)
			std_diff = np.mean(diff)
			print(f"Jupiter mass = {lab:5}Mj, Avrage deviation = {avg_diff:.5E}, std = {std_diff:.5E}")
		ax.tick_params(axis='both', which='major', labelsize=13)
		ax.set_xlabel("Time [yr]", fontsize=14)
		ax.set_ylabel("Deviation from circular orbit [AU]", fontsize=14)
		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16),
          ncol=4, fancybox=True, shadow=True, fontsize=14)
	plt.show()

