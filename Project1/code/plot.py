import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import platform
from scipy import stats

def run_program(methods, max_p):
	for method, p in zip(methods, max_p):
		# there are slight syntactic differences between Linux and python
		# concerning these commands; using platform-module to differentiate
		if platform.system() == 'Windows':
			cmd = ".\main.exe " + method + " " + str(p)
		else:
			cmd = "./main.exe " + method + " " + str(p)
		os.system(cmd)

def plot_data(methods, max_p):
	rel_error = []
	H = []
	for i, method in enumerate(methods):
		with sns.axes_style("darkgrid"):
			fig, ax = plt.subplots(nrows=1, ncols=1)
			ax.set(title=f"Calc vs analytical for {method} method", xlabel="x", ylabel="y")
		
		for j in range(1, max_p[i]+1):
			filename = "data/" + method + str(j) + ".txt"
			data = np.loadtxt(filename)
			x, calc, = data[:,0], data[:,1]

			h = 10**(-j)
			ax.plot(x, calc, label=f"h = {h}")

			if(method=="special"):
				eps = data[:,3]
				H.append(h)
				rel_error.append(max(eps))
		
		analytical = np.loadtxt(filename)[:,2]	
		ax.plot(x, analytical, label="Analytical", linestyle="dashed", c="k")
		ax.legend()
		plt.show()

	rel_error = np.array(rel_error)
	H = np.array(H)
	slope, const, r_value, p_value, std_err = stats.linregress(np.log10(H[:5]), np.log10(rel_error[:5]))
	print(slope, std_err)
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(yscale="log", xscale="log", xlabel="h", ylabel=r"Max($\epsilon$)")
		ax.plot(H[:5], 10**const * H[:5]**slope, c="k", \
				label="Linear fit", linestyle="dashed", marker='o', markersize=3)
		ax.scatter(H, rel_error, c="r", label="Computed")
		ax.legend()
		plt.show()

def plot_time(filename):
	filename = "data/" + filename + ".txt"
	data = np.loadtxt(filename)
	LU_max = 3
	thomas_max = 7
	LU = data[:,0:LU_max]/(1e6) # To mili seconds
	thomas = data[:,LU_max:(LU_max+thomas_max)]/(1e6) # To mili seconds
	thomas_sing = data[:,(LU_max+thomas_max):(LU_max+2*thomas_max)]/(1e6) # To mili seconds

	LU_middle = np.mean(LU, axis=0)
	thomas_middle = np.mean(thomas, axis=0)
	thomas_sing_middle = np.mean(thomas_sing, axis=0)

	LU_error = np.std(LU, axis=0)
	thomas_error = np.std(thomas, axis=0)
	thomas_sing_error = np.std(thomas_sing, axis=0)

	H_LU = [10**(-i) for i in range(1, LU_max+1)]	
	H_thomas = [10**(-i) for i in range(1, thomas_max+1)]

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(yscale="log", xscale="log", xlabel="h", ylabel=r"<t> [m/s]")
		ax.plot(H_LU, LU_middle, label="LU")
		ax.plot(H_thomas, thomas_middle, label="Thomas")
		ax.plot(H_thomas, thomas_sing_middle, label="Thomas singel valued")
		ax.legend()
		plt.show()

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(xscale="log", xlabel="h", ylabel=r"<t> [m/s]")
		ax.scatter(H_thomas, thomas_middle/thomas_sing_middle)
		plt.show()
if __name__ == "__main__":
	methods = ["LU", "Thomas", "Thomas_singval"]
	max_p = [1, 1, 7]
 
	# run_program(methods, max_p)
	# plot_data(methods, max_p)
	plot_time("times")
