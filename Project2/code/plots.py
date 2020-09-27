import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import sys
import platform
from scipy import stats

from file_reader import read_data_file

def run_bb():
	"""
	If there is no data file for BucklingBeam, run the code for N values
	"""
	N = [20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200]
	for n in N:
		os.system(".\BucklingBeam.exe " + str(n))



def plot_bb_eigvectors(run_index=0, vec_start=0, vec_end=0): 
	"""
	args:
		run_index: What run (from the BucklingBeam.dat file) to choose vectors from
				   preferably one with N > 100 (for better resolution in the plot)
		vec_start: The first eigenvector to plot
		vec_end: Up to and including this eigenvector to plot 
	"""
	runs = read_data_file("data/BucklingBeam.dat")
	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()

		if vec_start == vec_end: # If only one is needed
			vec_indexes = [0]
		else:
			vec_indexes = np.arange(vec_start, vec_end+1) # If multiple make a list with all desired Ns
		
		rho = np.linspace(0, 1, runs[run_index]('N')) 
		for i in vec_indexes:
			eigen_val = runs[run_index].vals[i]

			vecs = eigen_val*runs[run_index].vecs[:,i]
			
			ax.plot(rho, vecs, label=f"Eig vec: {i+1}")
		
		ax.set(xlabel=r"$\xi$", ylabel="HER SKAL DET STÅ NOE MEN VET IKKE HELT HVA")
		ax.legend()
		plt.show()



def plot_qo_groundstate(no_electrons, n): 
	"""
	args:
		no_electrons: either 'one' or 'two', specifying whether to plot data from
					  one- or two-electron systems.
		n			: energy level to plot (int)
	"""
	runs = read_data_file("data/QuantumOscillator_" + no_electrons + ".dat")
	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		
		for run in runs:
			N = run('N')
			rho_max = run('rho_max')

			lbl = f"N = {N}"
			lbl += r", $\rho_{max}$ = " + f"{rho_max}"
			if no_electrons == 'two':
				lbl += r", $\omega_r$ = " + f"{run('omega_r')}"
			lbl += ", $E_{%i}$ = " % n + f"{run.vals[n-1]}"

			ground_state = run.vecs[:,n-1] # Getting the eigenvector corresponding to the lowest eigenvalue
			ground_state *= np.sqrt((N+1)/rho_max) # Normalizing it
			rho_0 = rho_max / (N+1)
			rho = np.linspace(rho_0, rho_max - rho_0, N)

			ax.plot(rho, ground_state**2, label=lbl)

		ax.set(xlabel=r"$\rho$", ylabel="$u_{%i,0}(r)$" % n)
		ax.legend()
		plt.show()



def plot_convergence():
	"""
	Code to plot the convergence, also makes a linear model of the points
	"""
	runs = read_data_file("data/BucklingBeam.dat")

	N = []
	n_iter = [] 
	for run in runs:
		N.append(run('N')+1) # +1 since the the matrix is has N-1 x N-1 dims
		n_iter.append(run('n_iter')) # Number of iterations to reach epsilon


	slope, const, r_value, p_value, std_err = stats.linregress(np.log10(N), np.log10(n_iter))
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(xscale="log", yscale="log", xlabel="N", ylabel="Iterations before $\epsilon$")
		ax.plot(N, 10**const * N**slope, c="k" ,linestyle="dashed",\
				label=f"Linear fit, slope = {slope:.2f}$\pm${std_err:.2f}",marker='o', markersize=3)
		ax.scatter(N, n_iter)

		ax.set_xticks([20, 30, 40, 60, 100, 150, 200])
		ax.set_xticklabels(["$10$", "$20$", "$40$", "$60$", "$10^2$", r"$1.5 \times 10^2$", "$10^2$"])
		"""
		ax.set_yticks([7e2, 2.6e3, 1e4 ,3e4, 8.8e4])
		ax.set_yticklabels([r"$7 \times 10^2$", r"$2.6 \times 10^3$", "$10^4$", r"$3 \times 10^4$", r"$8.8 \times 10^4$"])
		"""
	ax.legend()
	plt.show()



def print_eigenvals(no_electrons, start_idx, stop_idx):
	"""
	This function prints a table of the eigenvalues for the specified data.
	args:
		no_electrons: either 'one' or 'two'; specifies whether to print the eigenvalues of the one- or two-electron systems
		start_idx   : (int) specifies the first eigenvalue to print, counting as 0, 1, ...
		stop_idx    : (int) specifies the last eigenvalue to print
	"""
	runs = read_data_file("data/QuantumOscillator_" + no_electrons + ".dat")

	no_cols = stop_idx - start_idx + 2
	rows = np.zeros((len(runs)+1, no_cols), dtype=object)

	rows[0,0] = r"Parameters [$N, \rho_\text{max}$]"
	rows[0,1:] = [r"$\lambda_" + f"{i}$" for i in np.arange(start=start_idx, stop=stop_idx+1)]
	for i, run in enumerate(runs):
		eigs = ["{:.3f}".format(eigval) for eigval in run.vals[start_idx:stop_idx+1]]
		rows[i+1,0] = f"${run('N')}$, ${run('rho_max')}$"
		rows[i+1,1:] = eigs

	print_table(rows)



def print_table(rows):
	"""
	This method prints a latex table of the elements given by the rows-matrix
	args:
		rows: (2D array-like of strings) The string elements that the table should include as (rows, cols).
	"""
	print(r"\begin{table}[H]")
	print(r"\begin{tabular}{" + "".join([r"|l" for i in range(rows.shape[1])]) + r"|}")
	print(r"\hline")

	for row in rows:
		row_txt = " & ".join(row)
		print("\t", row_txt, r" \\ \hline", sep="")

	print(r"\end{tabular}")
	print(r"\end{table}")



def get_flags():
	"""
	
	"""
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



def parse_flags(flags):
	"""
	Args:
		flags: Takes a list of of flags and parses them to choose what plots to show
	"""
	files = os.listdir("data/")

	if "h" in flags:
		print("The avalaible flags are")
		print("Not implemented yet...")
		sys.exit(1)

	if "v" in flags or "c" in flags: 
		if not "BucklingBeam.dat" in files:
			run_bb()
	if "v" in flags:
		plot_bb_eigvectors(4, vec_start=0, vec_end=1)

	if "c" in flags:
		plot_convergence()

	if "q" in flags:
		try:
			no_electrons = sys.argv[2]
			n = int(sys.argv[3])
		except IndexError:
			raise Exception("Using the quantum flag q requires the second argument to specify whether to plot \
							for 'one' electron or 'two' and that the third argument specifies the energy level \
							to plot for (1, 2, ...).")
		plot_qo_groundstate(no_electrons, n)

	if "e" in flags:
		try:
			no_electrons = sys.argv[2]
			start_idx = int(sys.argv[3])
			stop_idx = int(sys.argv[4])
		except IndexError:
			raise Exception("Using the quantum flag e requires the second argument to specify whether to plot \
							for 'one' electron or 'two' and that the third and fourth arguments specify first \
								and last eigenvalue to plot.")
		print_eigenvals(no_electrons, start_idx, stop_idx)



if __name__ == "__main__":
	
	flags = get_flags()
	
	check_compile()

	parse_flags(flags)