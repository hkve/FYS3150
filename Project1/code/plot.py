import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import platform
from scipy import stats

def run_program(methods, max_p):
	for method, p in zip(methods, max_p):
		# there are slight syntactic differences between Linux and Windows
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


if __name__ == "__main__":
	methods = ["LU", "Thomas", "special"]
	max_p = [1, 1, 7]

	# run_program(methods, max_p)
	plot_data(methods, max_p)

