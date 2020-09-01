import matplotlib.pyplot as plt
import numpy as np
import os as os

def run_program(methods, max_p):
	for i in range(len(methods)):
		cmd = "./main.exe " + methods[i] + " " + str(max_p[i]) 
		os.system(cmd)

def plot_data(methods, max_p):
	for i in range(len(methods)):
		for j in range(1, max_p[i]+1):
			fig, ax = plt.subplots(nrows=1, ncols=1)
			filename = "data/" + methods[i] + str(j) + ".txt"
			data = np.loadtxt(filename, delimiter=",")
			x, calc, analytical = data[:,0], data[:,1], data[:,2]

			h = 10**(-j)
			ax.set(title=f"Calc vs analytical for h = {h}", xlabel="x", ylabel="y")
			ax.plot(x, calc, label=methods[i])
			ax.plot(x, analytical, label="Analytical")
			plt.legend()
			plt.show()

"""
methods = ["sgeneral", "general"]
max_p = [3, 6]
"""
plot_data(methods, max_p)

