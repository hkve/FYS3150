import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os

def run_program(methods, max_p):
	for i in range(len(methods)):
		cmd = "./main.exe " + methods[i] + " " + str(max_p[i]) 
		os.system(cmd)

def plot_data(methods, max_p):
	rel_error = []
	H = []
	for i in range(len(methods)):
		with sns.axes_style("darkgrid"):
			fig, ax = plt.subplots(nrows=1, ncols=1)
			ax.set(title=f"Calc vs analytical for {methods[i]} method", xlabel="x", ylabel="y")
		
		for j in range(1, max_p[i]+1):
			filename = "data/" + methods[i] + str(j) + ".txt"
			data = np.loadtxt(filename)
			x, calc, = data[:,0], data[:,1]

			h = 10**(-j)
			ax.plot(x, calc, label=f"h = {h}")

			if(methods[i]=="special"):
				eps = data[:,3]
				H.append(h)
				rel_error.append(max(eps))
		
		analytical = np.loadtxt(filename)[:,2]	
		ax.plot(x, analytical, label="Analytical", linestyle="dashed", c="k")
		ax.legend()
		plt.show()

	rel_error = np.array(rel_error)
	H = np.array(H)
	print(H)
	fig, ax = plt.subplots(nrows=1, ncols=1)
	ax.set_yscale("log")
	ax.set_xscale("log")

	ax.scatter(H, rel_error)
	plt.show()
methods = ["sgeneral", "general", "special"]
max_p = [1, 1, 7]

plot_data(methods, max_p)

