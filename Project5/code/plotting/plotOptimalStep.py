def plotOptimalStep(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from optimalStep import findPercentAccepts, findOptimalStep
	from scipy.optimize import curve_fit
	from colour import Color
	import matplotlib as mpl
	import seaborn as sns
	omega = 1
	alphaStart = 0.1
	alphaEnd = 2.5
	dAlpha = 0.1

	stepStart = 0.1
	stepEnd = 5
	dStep = 0.01

	filename = "optimalStep1 copy.dat"

	if sim:
		os.system(f"../compiled/optimalStep.exe {5} {omega} {alphaStart} {alphaEnd} {dAlpha} {stepStart} {stepEnd} {dStep} {filename}")

	steps, alphas, percentAcceptAlphas = findPercentAccepts(filename)
	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(ncols=2, nrows=1,dpi = 125)
		ax = axes[0]
		ax2 = axes[1]
		ax.axhline(0.5, linestyle="dashed", c="k", zorder=10000, label=r"50% acceptance rate")
		ax.set_xlabel("Step length, $\delta$", fontsize=13)
		ax.set_ylabel("Acceptance rate", fontsize=13)
		ax.set_title(r"Acceptance rate for different $\alpha$")
		ax.set_xlim([0,5])


		colors = list(Color("blue", luminance=0.5).range_to(Color("red", luminance=0.45),len(alphas))) 
		colors = [color.get_rgb() for color in colors]
		
		for i,(alpha, percentAccepts) in enumerate(zip(alphas, percentAcceptAlphas)):	
			ax.plot(steps, percentAccepts,color=colors[i], alpha=0.6)
		cmap = mpl.colors.ListedColormap(colors)
		norm = mpl.colors.Normalize(vmin=alphas[0],vmax=alphas[-1])
		cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax2, orientation='vertical',ticks=[2.5, 2, 1.5, 1, 0.5] )
		cbar.ax.set_ylabel(r'$\alpha$', rotation=0, fontsize=18)
		ax.legend(loc='upper right', fancybox=True,shadow=True,fontsize=11)
		

		
		ax2.set_xlabel(rf"$\alpha$", fontsize=13)
		ax2.set_ylabel(r"Step size closest to 50% accepts", fontsize=13)
		optimalStep = findOptimalStep(steps, alphas, percentAcceptAlphas)
		for i,(alpha, oS) in enumerate(zip(alphas, optimalStep)):
			ax2.scatter(alpha, oS, color=colors[i])
		fitfunc = lambda x,s: s/x**(1/2)
		popt, pcov = curve_fit(fitfunc, alphas, optimalStep)
		print(popt[0])
		ax2.plot(alphas, fitfunc(alphas, *popt), color="k", alpha=0.7, lw=2,label=f"Fit: ${round(popt[0],2)}"+r"\cdot \alpha^{-1/2}$")
		ax2.legend(loc='upper right', fancybox=True,shadow=True,fontsize=11)
		ax2.set_title(r"Steps size yielding 50% acceptance rate")
		plt.show()

plotOptimalStep(sim = False)