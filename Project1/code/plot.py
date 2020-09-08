'''
This program acts as a main hub from which we can compile, run and plot all other scripts.
It also contains the methods for plotting various data written to file.

Executing this python file requires flags for the various operations.
-c : compiling the c++ programs
-r : running the methods that are called as arguments with c++
-t : running the time analysis
The flags can be combined, as in "-cr" to compile and run

To execute a method, and the method name (ie. LU) to the args, followed by the maximum
power of 10 for which to set n to as it solves the Poisson equation numerically.
To run multiple methods, you can just add arguments, alternating methods and max powers.
Run Example:
>>> python .\plot.py -cr LU 3

'''


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import sys
import platform
from scipy import stats



def plot_data(method, max_p):
	rel_error = []
	H = []
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(title=f"Calc vs analytical for {method} method", xlabel="x", ylabel="y")
	
	for j in range(1, max_p+1):
		filename = "data/" + method + str(j) + ".txt"
		data = np.loadtxt(filename)
		x, calc, = data[:,0], data[:,1]

		h = 10**(-j)
		ax.plot(x, calc, label=f"h = {h}")

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
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(yscale="log", xscale="log", xlabel="$h$", ylabel="Max($\epsilon$)")
		ax.plot(H[:5], 10**const * H[:5]**slope, c="k", \
				label=f"Linear fit, slope = {slope:.2f}$\pm${std_err:.2f}", linestyle="dashed",marker='o', markersize=3)
		ax.scatter(H, rel_error, c="r", label="Computed")
		ax.legend()
		plt.show()


def plot_time(filename):
	filename = "data/" + filename + ".txt"
	data = np.loadtxt(filename)
	LU_max = 3
	thomas_max = 7
	LU = data[:,0:LU_max]/(1e6) # To milliseconds
	thomas = data[:,LU_max:(LU_max+thomas_max)]/(1e6) # To milliseconds
	thomas_sing = data[:,(LU_max+thomas_max):(LU_max+2*thomas_max)]/(1e6) # To milliseconds

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
		ax.set(
			title = "Computation time of the algorithms",
			xscale = "log",
			yscale="log",
			xlabel = "$h$",
			ylabel = r"$\bar{t}$ [ms]"
		)
		ax.errorbar(H_LU, LU_middle, yerr=LU_error,label="LU")
		ax.errorbar(H_thomas, thomas_middle, yerr=thomas_error ,label="Thomas")
		ax.errorbar(H_thomas, thomas_sing_middle, yerr=thomas_sing_error, label="Thomas single valued")
		ax.legend()
		plt.show()

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(
			title = "Relative computation time",
			xscale = "log",
			xlabel = "$h$",
			ylabel = r"$\bar{t}_{8n} / \bar{t}_{4n}$"
		)
		thomas_diff = thomas_middle/thomas_sing_middle
		thomas_diff_error = thomas_diff*np.sqrt((thomas_error/thomas_middle)**2+(thomas_sing_error/thomas_sing_middle)**2) 

		for i in range(len(H_thomas)):
			ax.scatter(H_thomas[i], thomas_diff[i], c="r")
			ax.errorbar(H_thomas[i], thomas_diff[i], yerr=thomas_diff_error[i], c="r")
		plt.show()


		
if __name__ == "__main__":
	methods = ["LU", "Thomas", "Thomas_singval"] # available algorithms

	flags = [] # for the operations to execute
	try:
		if sys.argv[1][0] == "-": # check to see whether any flag is given
			for flag in sys.argv[1][1:]:
				flags.append(flag)
		else:
			raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")
	except IndexError:
		raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")


	if "h" in flags:
		print("The avalaible flags are:")
		print("-c : compile programs")
		print("-r : run programs")
		print("-t : plot time analysis")
		print("\nThe available methods are:")
		for method in methods:
			print(method)
		sys.exit(1)

		

	if "c" in flags: # compiles using makefile
		os.system("make compile")

	for method in methods:
		# checks to see if any argument matches a method
		if method in sys.argv:
			try:
				# tries to find the index of the current method in args
				# and set max_p to be the value after
				arg_idx = sys.argv.index(method)
				max_p = int(sys.argv[arg_idx+1])
			except ValueError:
				raise Exception("The argument following a method must be a valid int.")
			except IndexError:
				raise Exception("There was no argument following the method, a maximum p is required.")

			if "r" in flags: # runs the program
				os.system("make execute " + f"A={method} " + f"B={max_p}")

			# plots the results written to file
			plot_data(method, max_p)

	if "t" in flags: # plots the result of the time analysis
		plot_time("times")


