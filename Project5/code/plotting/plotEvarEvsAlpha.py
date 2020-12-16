def main(sim = False): 
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import seaborn as sns
	from colour import Color
	MCCs = 7
	alphaStart = 0.5
	alphaEnd = 1.5
	dAlpha = 0.01

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1
	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)
	colors = ["blue","red"][::-1]
	labels = ["With interactions", "Without interactions"]
	linestyles = ["solid", "dashed"]
	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=1, ncols=2, dpi=120)
		for i,mode in enumerate(["interactive", "noninteractive"]):
			# if i == 0:
			# 	continue
			if sim: #and i != 0:
				os.system(f"../compiled/variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} NoSave {mode}")

			E, E2 = np.loadtxt(f"../data/alpha{mode}.dat", unpack=True, usecols=(3,4))
			varE = E2-E*E			
			axes[0].plot(alphas, E,label=labels[i],color=colors[i], lw=4) #,label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
			Eargmin = np.argmin(E) 
			varEargmin = np.argmin(varE)
			varEargmin = list(alphas).index(0.88)

			print(alphas[varEargmin],varE[varEargmin])
			
			axes[1].plot(alphas, varE,label=labels[i],color=colors[i], lw=3.5, alpha=0.8) #,label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
			axes[0].scatter(alphas[Eargmin], E[Eargmin],  label=f"$\\alpha = {alphas[Eargmin]}, E = {round(E[Eargmin],2)}$", color=Color(colors[i],luminance=0.33).get_rgb(), zorder=100)

		for i in range(2):
			axes[i].set_xlabel(r"$\alpha$", fontsize=16)
		axes[0].set_ylabel(r"$E_{L1}$ [a.u]", fontsize=16)
		axes[1].set_ylabel(r"$\sigma_{E_{L1}}^2$ [a.u$^2$]", fontsize=18)
		axes[0].set_title(r"Energy as function of $\alpha$", fontsize=17)
		axes[1].set_title(r"Variance energy as function of $\alpha$. $\omega = 1$.", fontsize=17)
		lines, labels = axes[0].get_legend_handles_labels()
		fig.legend(lines, labels,loc='upper center', fancybox=True, ncol=2,shadow=True,fontsize=14)
		#axes[0].legend()
		plt.show()


