def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import seaborn as sns
	from subprocess import run

	Tstart = 2
	Tend = 2.5
	dT = 0.05
	MCCs = int(1e5)
	#L_ = [40, 60, 80, 100]
	L_ = [5, 10, 15, 20]
	filenames = [f"Exp_values_L{L}.dat" for L in L_]

	if sim:
		#os.system("../cpp/Paralell.out 5 1000000 2 3 0.1 data.out")
		for L, filename in zip(L_, filenames):
			if filename in os.listdir("../data/"):
				os.system(f"rm ../data/{filename}")

			os.system(f"../cpp/paralell.out {L} {MCCs} {Tstart} {Tend} {dT} {filename}")

	fig, axes = plt.subplots(nrows=1, ncols=2)
	with sns.axes_style("darkgrid"):
		for L, filename in zip(L_, filenames):
			data = np.loadtxt(f"../data/{filename}")
			varE, varM, T, = data[:,5], data[:,6], data[:,7]

			Cv = varE/T**2
			X = varM/T
			axes[0].plot(T, Cv, label="L = {L}")
			axes[1].plot(X, T, label="L = {}")

		axes[1].legend()
		axes[2].legend()

	plt.show()

main(True)