def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from scipy.optimize import curve_fit
	import seaborn as sns
	modes = ["noninteractive", "interactive"]
	filenames = [f"virial_{mode}.dat" for mode in modes]

	MCCs = 7
	omegaStart = 0.01 
	omegaEnd = 1
	N_omega = 20

	if sim:
		for mode, filename in zip(modes, filenames):
			#if filename in os.listdir("../data/"):
			#	os.system(f"rm ../data/{filename}")
			#os.system(f"../compiled/virialTheorem.exe {MCCs} {omegaStart} {omegaEnd} {N_omega} {mode}")
			pass
	colors = ["red", "blue"][::-1]
	labels = ["Without interactions", "With interactions"]
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(dpi=140)
		#ax.set_ylim(0,4)

		for i,(mode, filename) in enumerate(zip(modes, filenames)):
			omega, E, EE, r12, V, r12_inverse = np.loadtxt("../data/virial_noninteractive.dat", unpack=True)
		
			K = E-V
			if mode == "interactive":
				V += r12_inverse
			
			ax.scatter(omega, V/K, label=labels[i], color=colors[i], zorder=100, alpha=0.8)

			if i == 1:
				fitfunc = lambda x,m: 1+m/x**(0.5)
				w_hires = np.linspace(omega[0],omega[-1], 10000, endpoint=True)
				popt, pcov = curve_fit(fitfunc, omega, V/K)
				print(popt)
				ax.plot(w_hires, fitfunc(w_hires, *popt), color="k", alpha=0.7, lw=2, label=f"Fit: $1+ 0.532\\cdot \\omega"+r"^{-0.5}$", linestyle="--")


		
		ax.set_xlabel(r"$\omega$", fontsize=14)
		ax.set_ylabel(r"$\langle V \rangle / \ \langle T \rangle$", fontsize=14)
		ax.legend(fancybox=True, shadow=True, loc="upper right", fontsize=12)
		ax.set_title(r"Ratio $\langle V \rangle / \ \langle T \rangle$ for both systems", fontsize=14)
		plt.show()
	
