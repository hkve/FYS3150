import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file

def plot_Sun_Earth_Jupiter(dt=0.001, T = 10, jupiter_scale = 1):
	body_dict = {
				"Sun": [0,0,0,0,0,0],
				"Earth": [1,0,0,0,2*np.pi,0],
				"Jupiter": [-5.2,0,0,0,-2*np.pi/np.sqrt(5.2),0]
	}

	initFilename = f"SunEarthJupiter_stable_{jupiter_scale}m.dat"
	outFilename = f"SunEarthJupiter_stable_{jupiter_scale}m.dat"
	setInitialConditions(initFilename, body_dict, fixedCoM=False)
	
	N = int(T/dt)
	N_write = 10000
	if True:
		master_call = f"python3 master.py -method verlet -sys initData/{initFilename} \
							-out {outFilename} -Nwrite {N_write} -time years --fixSun {dt} {N}" 
		subprocess.call(master_call.split())

	system = read_data_file(outFilename)
	
	
	rE, rJ = system["Earth"].r, system["Jupiter"].r
	
	fig, ax = plt.subplots()
	ax.plot(rE[0], rE[1])
	ax.plot(rJ[0], rJ[1])
	plt.show()
	
plot_Sun_Earth_Jupiter()