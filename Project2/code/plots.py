import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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



def run_qo2(slash, N, rho_max_list):
	"""
	Run the QuantumOscillator function for two electrons for the relevant omega_rs 
	and the given parameters N and rho_max.s
	"""
	open("data/QuantumOscillator_two.dat", "w").close() # empty the contents of the .dat file
	omega_r_list = [0.01, 0.5, 1, 5]

	for omega_r, rho_max in zip(omega_r_list, rho_max_list):
		args = " ".join(["2", str(N), str(rho_max), str(omega_r)])
		os.system(slash + "QuantumOscillator.exe " + args)


def plot_error_vs_rho_N(slash):
	"""
	Code to produce 3D plot of error vs rho_max and N
	Demonstrates how the error (in this case 4 lowest eigenvalues)
	is dependent on both rho_max and N
	"""
	filename = "data/QuantumOscillator_one.dat"
	
	if filename.replace("data/", "") in os.listdir("data/"): # If a file is already there, delete it
		os.remove(filename)
	
	N_start, N_end = 100, 200 				# start and stop for the matrix sices
	rho_max_start, rho_max_end = 3.5, 6.5   # start and stop for the rho_maxes to test
	dN = 5 									# space between each matrix sice
	dRho_max = 0.1							# space between each rho_max

	N = np.arange(N_start, N_end+dN, dN) 								# array of all Ns
	Rho_max = np.arange(rho_max_start, rho_max_end+dRho_max, dRho_max)  # array of all rho_maxes
	error = np.zeros((len(N), len(Rho_max)))							# grid to hold error
	
	for n in N: # Run the c++ code
		for rho_max in Rho_max:
			args = "1 " + str(n) + " " + str(rho_max) + " 1"
			os.system(slash + "QuantumOscillator.exe " + args)
	
	
	runs = read_data_file(filename, drop_vecs=True) # Read .dat file, excluding eigenvectors due to memory
	ana_val = np.array([3, 7, 11, 15])              # 4 lowest analytical eigenvalues
	
	i, j = 0, 0 # i is the index for each N, j is index for each rho_max
	for run in runs:
		num_val = run.vals[:4] # Get 4 lowest eigenvalues
		err = np.max(np.abs(num_val-ana_val)) # Calc max error of the 4 eigenvalues
		error[i,j] = err # Store in error grid

		j += 1 # Go to next rhomax for a N 
		if j == len(Rho_max): # If we have done all rhomax for this N
			j = 0
			i += 1 # Go to next N
	
	Rho_max, N = np.meshgrid(Rho_max, N) 
	error = np.log10(error) # Log10 of error due to visibility

	fig = plt.figure()
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(N, Rho_max, error,cmap="gnuplot")
	ax.set_xlabel("N", fontsize=13)
	ax.set_ylabel(r"$\rho_{max}$", fontsize=13)
	ax.set_zlabel("log10($|\lambda_{num}-\lambda_{ana}|$)", fontsize=13)
	plt.show()



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



def plot_qo_eigvecs(no_electrons, n): 
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
			lbl += ", $\lambda_{%i}$ = " % n + f"{run.vals[n-1]:.2f}"

			ground_state = run.vecs[:,n-1] # Getting the eigenvector corresponding to the lowest eigenvalue
			ground_state *= np.sqrt((N+1)/rho_max) # Normalizing it
			rho_0 = rho_max / (N+1)
			rho = np.linspace(rho_0, rho_max - rho_0, N)

			ax.plot(rho, ground_state**2, label=lbl)

		ax.set_xlabel(r"$\rho$", fontsize=12)
		ax.set_ylabel("$|u_{%i,0}(r)|^2$" % n, fontsize=12)
		ax.legend()
		plt.show()



def plot_time_difference(slash):
	"""
	Plotting the relationship of Jacobi time and Armadillo time vs matrix sizes (N)
	y-axis is logarithmic to show more details for low N
	"""
	filename = "data/time.dat"
	N = np.array([i*2 for i in range(1,101)]) # Testing Ns from 2 to 200 with dN = 2

	if  not filename.replace("data/", "") in os.listdir("data/"): # if the file is not there
		for n in N:
			os.system(slash + "BucklingBeam.exe " + str(n) + " 1") # "1" to indicate that we want to run timing

	data = np.loadtxt(filename) # Load from .dat file
	time_jacobi, time_arma = data[:,0], data[:,1] 
	relative = time_jacobi/time_arma

	line_i = np.argmin(np.abs(relative-1)) # The index where the relationship is closest to 1 
										   # I.e. the matrix size where the two methods are approx equal
	
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
		eigs = ["${:.4f}$".format(eigval) for eigval in run.vals[start_idx:stop_idx+1]]
		rows[i+1,0] = f"${run('N')}$, ${run('rho_max')}$"
		rows[i+1,1:] = eigs

	print_table(rows)



def print_table(rows):
	"""
	This method prints a LaTeX table of the elements given by the rows-matrix
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
	if platform.system() == "Windows":
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

	if "v" in flags:
		plot_bb_eigvectors(4, vec_start=0, vec_end=1)

	if "c" in flags:
		plot_convergence()

	if "t" in flags: 
		plot_time_difference(slash)
	
	if "r" in flags:
		plot_error_vs_rho_N(slash)

	if "q" in flags:
		try:
			no_electrons = sys.argv[2]
			n = int(sys.argv[3])

		except IndexError:
			raise Exception("Using the quantum flag q requires args as: type of system ('one' or 'two'), eigvec-level (1,2,...)")
		if no_electrons == 'two':
			try:
				N = int(sys.argv[4])
				rho_max_list = [float(arg) for arg in sys.argv[5:9]]
			except IndexError:
				raise Exception("For the two-electron system further args are required as:\n" + \
								"N (number of integration points), rho_maxes (4 rho_maxes to run with omega_r = 0.01, 0.5, 1, 5)")
			run_qo2(slash, N, rho_max_list)
		plot_qo_eigvecs(no_electrons, n)

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
	
	check_compile()
	
	flags = get_flags()

	parse_flags(flags)