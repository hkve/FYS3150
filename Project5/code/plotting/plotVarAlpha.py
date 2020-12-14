import matplotlib as mpl
from colour import Color

def varAlpha(sim = False): 
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import seaborn as sns

	MCCs = 7
	alphaStart = 0.5
	alphaEnd = 1.5
	dAlpha = 0.05

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1
	

	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)
	colors = list(Color("pink", luminance=0.5).range_to(Color("green", luminance=0.5),N_alpha)) 
	colors = [color.get_rgb() for color in colors]
	for filename, mode,title in zip(["varAlphaNoninteractive", "varAlphaInteractive"], ["noninteractive", "interactive"], ["Without interactions", "With interactions"]):
		if sim:
			os.system(f"../compiled/variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} Save {mode}")
		with sns.axes_style("darkgrid"):
			fig, axes = plt.subplots(nrows=1, ncols=2, dpi=120)
			axes[0].set(xscale="log")
			axes[1].set(xscale="log")
			for i, alpha in enumerate(alphas):
				data = np.loadtxt(f"../data/{filename}{i}.dat")
				cyles, E, E2 = data[:,0], data[:,1], data[:,2]

				varE = E2-E*E
				zorder = 1#int(10-10*np.abs(alpha -1))
				axes[0].plot(cyles, E, color=colors[i], zorder=zorder, lw=1.7,alpha=0.72)#, label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
				axes[1].plot(cyles, varE, color=colors[i], zorder=zorder,lw=1.7,alpha=0.72)#, label=fr"$\alpha$ = {alpha:.1f} $\sigma_E^2$ = {varE[-1]:.2f}")
			cmap = mpl.colors.ListedColormap(colors)
			norm = mpl.colors.Normalize(vmin=alphas[0],vmax=alphas[-1])
			cbar = axes[1].figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[1], orientation='vertical',ticks=[1.5, 1.25, 1, 0.75, 0.5] )

			cbar.ax.set_ylabel(r'  $\alpha$', rotation=0, fontsize=14)
			axes[0].set_title(f"Energy, {title}")
			axes[1].set_title(f"Variance of energy, {title}")

			axes[0].set_ylabel(r"$E_{L1}$ $[\omega]$", fontsize=15)
			axes[0].set_xlabel(r"$\tau$", fontsize=15)
			axes[1].set_ylabel(r"$\sigma_{E_{L1}}^2$  $[\omega^2]$", fontsize=14)
			axes[1].set_xlabel(r"$\tau$", fontsize=15)
			plt.show()
			


def varAlphaInteractive(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	MCCs = 7
	alphaStart = 0.85
	alphaEnd = 1.15
	dAlpha = 0.05

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1

	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)

	if sim:
		os.system(f"../compiled/variateAlpha.exe {5} {alphaStart} {alphaEnd} {dAlpha} Save interactive")

	fig, axes = plt.subplots(nrows=1, ncols=2)
	axes[0].set(xscale="log")
	axes[1].set(xscale="log")
	for i, alpha in enumerate(alphas):
		data = np.loadtxt(f"../data/varAlphaInteractive{i}.dat")
		cyles, E, E2 = data[:,0], data[:,1], data[:,2]

		varE = E2-E*E
		skip = 0
		axes[0].plot(cyles[skip:], E[skip:], label=fr"$\alpha$ = {alpha:.2f}, E = {E[-1]:.2f}")
		axes[1].plot(cyles[skip:], varE[skip:], label=fr"$\alpha$ = {alpha:.2f} $\sigma_E^2$ = {varE[-1]:.2f}")

	axes[0].legend()
	axes[1].legend()
	plt.show()

varAlpha(sim=True)
#varAlphaInteractive(sim = True)