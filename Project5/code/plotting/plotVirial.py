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
			omega, E, EE, r12, V, r12_inverse = np.loadtxt(f"../data/{filename}", unpack=True)
			print(np.mean(E), filename)
		
			if mode == "interactive":
				V += r12_inverse

			K = E-V
			ax.scatter(omega, K/V, label=labels[i], color=colors[i], zorder=100, alpha=0.8)

			if i == 1:
				alpha = 0.994
				fitfunc = lambda x,b: 1+b/(alpha*x)**(0.25) #passet best for 1/r12. Gir b = 0.6
				#fitfunc = lambda x,b: 2*np.log(2)+b/(alpha*x)**(0.5) # passet best for r12_inverse. Gir b = 0.47024 = ln(1+0.6)
				#fitfunc = lambda x,a,b,c: a+b/(alpha*x)**(c) 
				fitfunc = lambda x,b,c: 1+b/(alpha*x)**(c)

				popt, pcov = curve_fit(fitfunc, omega, K/V)
				b,c = tuple(popt)
				print(f"Zero point: {(-b)**(1/c)}")
				#ax.plot(w_hires, fitfunc(w_hires, *popt), color="k", alpha=0.7, lw=2, label=f"Fit: $1+ 0.532\\cdot \\omega"+r"^{-0.5}$", linestyle="--")


		
		ax.set_xlabel(r"$\omega$ [a.u]", fontsize=14)
		ax.set_ylabel(r"$\langle T \rangle / \ \langle V \rangle$", fontsize=14)
		ax.legend(fancybox=True, shadow=True, loc="lower right", fontsize=12)
		ax.set_title(r"Ratio $\langle T \rangle / \ \langle V \rangle$ for both systems", fontsize=14)
		plt.show()
	
