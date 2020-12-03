def varAlpha(sim = False): 
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import seaborn as sns

	#MCCs = 7
	alphaStart = 0.5
	alphaEnd = 1.5
	dAlpha = 0.01

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1
	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)
	colors = ["blue","red"][::-1]
	r12_ = []
	labels = ["With interactions", "Without interactions"]
	linestyles = ["solid", "dashed"]
	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=1, ncols=2, dpi=120)
		for i,mode in enumerate(["interactive", "noninteractive"]):
			if sim and i != 0:
				os.system(f"../compiled/variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} {mode} NoSave")

			E, r12 = np.loadtxt(f"../data/alpha{mode}.dat", unpack=True, usecols=(3,5))
			r12_.append(r12)
			
			axes[0].plot(alphas, E,label=labels[i],color=colors[i], lw=4) #,label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
			axes[1].plot(alphas, r12,label=labels[i],color=colors[i], lw=2+2.2-1.5*i,alpha=0.4*i/2+0.6, linestyle=linestyles[i]) #,label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")

		for i in range(2):
			axes[i].set_xlabel(r"$\alpha$", fontsize=16)
		axes[0].set_ylabel(r"$E_{L1}$ $[\omega]$", fontsize=16)
		axes[1].set_ylabel(r"$r_{12}$", fontsize=18)
		axes[0].set_title(r"Energy as function of $\alpha$", fontsize=17)
		axes[1].set_title(r"Distance $r_{12}$ as function of $\alpha$", fontsize=17)
		lines, labels = axes[0].get_legend_handles_labels()
		fig.legend(lines, labels,loc='upper center', fancybox=True, ncol=2,shadow=True,fontsize=14)
		#axes[0].legend()
		plt.show()



varAlpha(sim = False)