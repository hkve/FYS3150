import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import sys
import platform
from scipy import stats

from file_reader import read_data_file

def run_bb(slash):
	"""
	If there is no data file for BucklingBeam, run the code for N values
	"""
	N = [20, 25, 30, 35, 40, 50, 60, 70, 80, 100, 120, 140, 160, 180, 200]
	for n in N:
		os.system(slash + "BucklingBeam.exe " + str(n))

def run_qo1(slash):
	"""
	If there is no data file for QuantumOscillator_one, run the code for rho_max and N values
	"""
	N = [100, 150, 200]
	rho_max = [4 + 0.1*i for i in range(21)]

	for rho in rho_max:
		args = "1 " + str(N[0]) + " " + str(rho) + " 1"
		os.system(slash + "QuantumOscillator.exe " + args)

def plot_bb_eigvectors(run_index=0, vec_start=0, vec_end=0): 
	"""
	args:
		run_index: What run (from the BucklingBeam.dat file) to choose vectors from
				   preferably one with N > 100 (for better resolution in the plot)
		vec_start: The first eigenvector to plot
		vec_end  : Up to and including this eigenvector to plot 
	"""
	runs = read_data_file("data/BucklingBeam.dat")
	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()

		if vec_start == vec_end: # If only one is needed
			vec_indexes = [0]
		else:
			vec_indexes = np.arange(vec_start, vec_end+1) # If multiple make a list with all desired Ns
		
		N = runs[run_index]("N")
		rho = np.linspace(0, 1, N) 
		ana_vec = np.zeros((N, vec_end-vec_start+1)) # Array for the analytical eigenvectors
		
		for j in range(vec_end-vec_start+1): # For each vector
			for i in range(N): # For each vector element
				ana_vec[i,j] = np.sin((j+1)*(i+1)*np.pi/(N+1))
				
			ana_vec[:,j] /= np.linalg.norm(ana_vec[:,j]) # Normalise the j'th eigenvector
		
		for i in vec_indexes: # For each eigenvector
			num_vecs = runs[run_index].vecs[:,i] # Numerical eigenvectors
			eigen_val = runs[run_index].vals[i]  # Analytical eigenvectors

			
			ax.plot(rho, num_vecs, label=f"Eigen Vec: {i+1}")
			ax.plot(rho, ana_vec[:,i], c="k", \
					linestyle="dashed", dashes=(5,10))

		ax.set_xlabel(xlabel=r"$\xi$", fontsize=14)
		ax.set_ylabel(ylabel=r"$u(\xi)$", fontsize=14)
		ax.legend(fontsize=12)
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



def plot_rho_max(no_electrons, n):
	runs = read_data_file("data/QuantumOscillator_" + no_electrons + ".dat")

	analytical_eig = np.array([3, 7, 11]) # Analytical result for 3 lowest eigenvalues
	n_runs = len(runs)
	max_error = np.zeros(n_runs)
	rho_max = np.zeros(n_runs)

	for i in range(n_runs):
		error = abs(runs[i].vals[0:3]-analytical_eig)
		max_error[i] = np.max(error)	
		rho_max[i] = runs[i]("rho_max")

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(xlabel=r"$\rho_{max}$", ylabel="$|\lambda_{num}-\lambda_{ana}|$")
		ax.plot(rho_max, max_error, c="k", linestyle="dashed")

		ideal_rho_i = np.argmin(max_error)
		ideal_rho = abs(rho_max[ideal_rho_i])
		ideal_error = abs(max_error[ideal_rho_i])
		ideal_rho_error = abs(ideal_rho-rho_max[ideal_rho_i+1])
		plt.scatter(ideal_rho, ideal_error, c="r", s=100, \
					label=r"$\rho_{ideal} = $" + f"{ideal_rho} $\pm$ {ideal_rho_error:.1f}")
		
	plt.legend()
	plt.show()

def plot_time_difference(slash):
	filename = "data/time.dat"
	N = np.array([i*2 for i in range(1,101)])

	if  not filename.replace("data/", "") in os.listdir("data/"):
		for n in N:
			os.system(slash + "BucklingBeam.exe " + str(n) + " 1")

	data = np.loadtxt(filename)
	time_jacobi, time_arma = data[:,0]/1e9, data[:,1]/1e9
	relative = time_jacobi/time_arma

	line_i = np.argmin(np.abs(relative-1))
	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(yscale="log")
		ax.set_xlabel("N", fontsize=17)
		ax.set_ylabel(r"$t_{i}/t_{a}$", fontsize=17)
		ax.hlines(relative[line_i], N[0], N[-1], linestyle="dashed", label=r"$t_{j} \approx t_{a}$")
		ax.plot(N, relative)
		ax.legend(fontsize=14)
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
		ax.set(xscale="log", yscale="log")
		ax.set_xlabel("N", fontsize=14)
		ax.set_ylabel("Iterations before $\epsilon$", fontsize=14)
		ax.plot(N, 10**const * N**slope, c="k" ,linestyle="dashed",\
				label=f"Linear fit, slope = {slope:.2f}$\pm${std_err:.2f}",marker='o', markersize=3)
		ax.scatter(N, n_iter)

		ax.set_xticks([20, 30, 40, 60, 100, 150, 200])
		ax.set_xticklabels(["$10$", "$20$", "$40$", "$60$", "$10^2$", r"$1.5 \times 10^2$", "$10^2$"])
		"""
		ax.set_yticks([7e2, 2.6e3, 1e4 ,3e4, 8.8e4])
		ax.set_yticklabels([r"$7 \times 10^2$", r"$2.6 \times 10^3$", "$10^4$", r"$3 \times 10^4$", r"$8.8 \times 10^4$"])
		"""
		ax.legend(fontsize=12)
	
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
	rows[0,1:] = [r"$\lambda_" + f"{i}$" for i in np.arange(start=start_idx+1, stop=stop_idx+2)]
	for i, run in enumerate(runs):
		eigs = ["{:.3f}".format(eigval) for eigval in run.vals[start_idx:stop_idx+1]]
		rows[i+1,0] = f"${run('N')}$, ${run('rho_max')}$"
		rows[i+1,1:] = eigs

	print_table(rows)



def print_table(rows):
	"""
	This method prints a LaTex table of the elements given by the rows-matrix
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
	Reads off the flags given for execution.
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

	program_names = ["main.exe", "BucklingBeam.exe", "QuantumOscillator.exe"] # Program names
	program_makefile_names = ["compile_main", "compile_bb", "compile_qo"] # Program names in the makefile

	for i, name in enumerate(program_names): # Loops over all programs that should be compiled
		if not name in files: # If missing, compile it
			os.system("make " + program_makefile_names[i]) 
			print(f"Compiled {program_names[i]}")



def parse_flags(flags):
	"""
	Args:
		flags: Takes a list of of flags and parses them to choose what plots to show
	"""
	slash = "./"
	if platform.platform() == "Windows":
		slash = ".\\"

	files = os.listdir("data/")

	if "h" in flags:
		print("The avalaible flags are")
		print("-v\t Plot two eigenvectors Buckling Beam")
		print("-c\t Plot number of iteration as a function of N")
		print("-r\t Plot error against different rho values for one electron")
		sys.exit(1)

	# Making data files
	if "v" in flags or "c" in flags: 
		if not "BucklingBeam.dat" in files:
			run_bb(slash)
	if "r" in flags:
		if not "QuantumOscillator_one.dat" in files:
			run_qo1(slash)

	if "v" in flags:
		plot_bb_eigvectors(4, vec_start=0, vec_end=1)

	if "c" in flags:
		plot_convergence()

	if "r" in flags:
		plot_rho_max("one", 10)

	if "q" in flags:
		try:
			no_electrons = sys.argv[2]
			n = int(sys.argv[3])
		except IndexError:
			raise Exception("Using the quantum flag q requires the second argument to specify whether to plot "\
							+ "for 'one' electron or 'two' and that the third argument specifies the energy level" \
							" to plot for (1, 2, ...).")
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

	if "t" in flags: 
		plot_time_difference(slash)


if __name__ == "__main__":
	
	check_compile()
	
	flags = get_flags()

	parse_flags(flags)