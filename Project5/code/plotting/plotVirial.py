def plotVirial(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	modes = ["noninteractive", "interactive"]
	filenames = [f"virial_{mode}.dat" for mode in modes]

	MCCs = 7
	omegaStart = 0.01
	omegaEnd = 1
	N_omega = 20

	if sim:
		for mode, filename in zip(modes, filenames):
			if filename in os.listdir("../data/"):
				os.system(f"rm ../data/{filename}")
			os.system(f"../compiled/virialTheorem.exe {MCCs} {omegaStart} {omegaEnd} {N_omega} {mode}")
	
	fig, ax = plt.subplots()
	#ax.set_ylim(0,4)
	for mode, filename in zip(modes, filenames):
		omega, E, EE, r12, V, r12_inverse = np.loadtxt("../data/virial_noninteractive.dat", unpack=True)
	
		K = E-V
		if mode == "interactive":
			V += r12_inverse
		

		ax.scatter(omega, V/K, label=f"{mode}")
	
	ax.set(xlabel=r"$\omega$", ylabel=r"$\langle V \rangle / \ \langle T \rangle$")
	ax.legend()
	plt.show()
	
plotVirial(sim=True)