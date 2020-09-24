import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import sys
import platform
from scipy import stats

from file_reader import read_data_file

def run_bb():
	N = [(10+20*i) for i in range(0,25)]
	for n in N:
		os.system("./BucklingBeam.exe " + str(n))

def plot_error():
	runs = read_data_file("BucklingBeam.txt")

	# calculate analytical eigenvalues
	vals = []
	H = []
	max_error = []
	for i in range(len(runs)):
		N = runs[i].N
		vals = np.zeros(N, dtype=float)
		h = 1/N
		d = 2/h**2
		a = -1/h**2

		H.append(h)
		for j in range(1,N+1):
			vals[j-1] = d+2*a*np.cos(j*np.pi/N)
		
		rel_error = np.abs((runs[i].vals-vals)/vals)
		max_error.append(np.max(rel_error))
		
	slope, const, r_value, p_value, std_err = stats.linregress(np.log10(H), np.log10(max_error))
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=1)
		ax.set(yscale="log", xscale="log", xlabel="$h$", ylabel="Max($|\epsilon_r$)")
		ax.plot(H, 10**const * H**slope, c="k", \
				label=f"Linear fit, slope = {slope:.3f}$\pm${std_err:.3f}", linestyle="dashed",marker='o', markersize=3)
		ax.scatter(H, max_error, c="r", label="Computed")
		ax.legend()
		plt.show()


def get_flags():
	flags = []
	try:
		if sys.argv[1][0] == "-": # check to see whether any flag is given
			for flag in sys.argv[:][1:]:
				flags.append(flag[:][1:])
		else:
			raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")
	except IndexError:
		raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")

	return flags

def check_compile():
	"""
	Function for compiling needed programs in case of missing .exe files
	"""
	files = os.listdir() # Gets filenames from directory 

	if not "data" in files: # Check for missing directory to store datafiles
		os.mkdir("data")
		print("Created missing data directory")

	program_names = ["main.exe", "BucklingBeam.exe"] # Program names
	program_makefile_names = ["compile_main", "compile_bb"] # Program names in the makefile

	for i, name in enumerate(program_names): # Loops over all programs that should be compiled
		if not name in files: # If missing, compile it
			os.system("make " + program_makefile_names[i]) 
			print(f"Compiled {program_names[i]}")

def prase_flags(flags):
	os.chdir("data/")
	files = os.listdir()

	if "h" in flags:
		print("The avalaible flags are")
		print("Not implemented yet...")
		sys.exit(1)

	if "e" in flags:
		if not "BucklingBeam.txt" in files:
			run_bb()
		else:
			plot_error()

if __name__ == "__main__":
	
	flags = get_flags()
	
	check_compile()

	prase_flags(flags)