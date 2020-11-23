def main(sim=False):
	import numpy as np
	import matplotlib.pyplot as plt
	import seaborn as sns
	import os

	Lstart = 20
	Lend = 100
	dL = 5
	MCCs = 5
	Nthreads = 4
	Ntests = 3*Nthreads
	methods = ["Series", "Parallel"]

	filenames = [f"timeComparison_{Lstart}to{Lend}_dL{dL}_N{Ntests}_MCC{MCCs}_{method}.dat" for method in methods]

	MCCs = int(10**MCCs)

	if sim:
		for method, filename in zip(methods, filenames):
			os.chdir("../cpp/")
			if filename in os.listdir("../data/"):
				print(filename)
				os.system(f"rm ../cpp/{filename}")
			if method == "Series":
				os.system("g++ -O3 -o timeComparison.out timeComparison.cpp IsingModel2DParalell.cpp")
			if method == "Parallel":
				os.system("g++ -O3 -o timeComparison.out timeComparison.cpp IsingModel2DParalell.cpp -fopenmp")

			os.system(f"./timeComparison.out {Lstart} {Lend} {dL} {MCCs} {Ntests} {Nthreads} {filename}")
			print(f"Done L=[{Lstart},{Lend}] for {method}, {Ntests} tests for each L")
	os.chdir("../python/")

	SeriesTime = []
	ParallelTime = []

	with sns.axes_style("darkgrid"):
		fig, ax  = plt.subplots()
		ax.set_xlabel("L", fontsize=14)
		ax.set_ylabel("Time used [s]", fontsize=14)
		ax.tick_params(axis='both', which='major', labelsize=13)
		colors = ["teal", "darkmagenta"]
		for method, filename, color in zip(methods, filenames, colors):
			data = np.loadtxt(f"../data/{filename}")
			L, time = data[:,0], data[:,1]

			ax.plot(L, time, label=f"{method}", marker="o", linestyle="dashed", color=color, linewidth=2)

			if method == "Series":
				SeriesTime = time
			if method == "Parallel":
				ParallelTime = time

		ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.16),
	          ncol=2, fancybox=True, shadow=True, fontsize=15)

		for L_, series, parallel in zip(L, SeriesTime, ParallelTime):
			ratio = series/parallel
			print(f"L = {L_:.0f}, Series = {series:.4f}s, Parallel={parallel:.4f}s, ratio = {ratio:.4f}")

		plt.show()

main(sim=False)