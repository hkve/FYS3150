import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

from getInitialConditions import setInitialConditions
from file_reader import read_data_file

def plot_circular_orbit(dt, T_end, method):
	"""
	Makes plot of stable Sun/Earth orbit
	Sun at origin no init vel, earth x = 1 AU vx = 2pi * AU/yr

	Args:
		dt: (float) time step in days
		T_end: (int/float) end time in years
		method: "euler" or "verlet"
	"""
	d2y = 365  	# Day to year conversion
	N = int(T_end*d2y/dt)

	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,2*np.pi/d2y,0,0]}
	
	initFilename = "SunEarthStable_init.dat"
	outFilename = "SunEarthStable_" + str(T_end) + "_" + str(N)  + ".dat"

	initDataDIR = os.listdir("initData")
	if not initFilename in initDataDIR: 				# Check if init file exists
		setInitialConditions(initFilename, body_dict)
	initFilename = "initData/" + initFilename 
	
	dataDIR = os.listdir("data")
	if not outFilename in dataDIR: 						# Check if data file exists
		master_call = f"python3 master.py -method {method} -sys {initFilename} \
						-out {outFilename} {dt} {N}" 
		subprocess.call(master_call.split())

	
	system = read_data_file(outFilename)

	r = system["Earth"].r
	plt.plot(r[0],r[1])
	plt.show()
	
plot_circular_orbit(dt=0.1, T_end=1, method="euler")