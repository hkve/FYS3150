def main(sim = False, omega=1,a0=0.9, a1=3, da=0.05, step0=0.5,step1=3, dstep=0.01, MCCs=6):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from optimalStep import findPercentAccepts, findOptimalStep
	from scipy.optimize import curve_fit
	from colour import Color
	import matplotlib as mpl
	import seaborn as sns
	alphaStart = a0
	alphaEnd = a1
	dAlpha = da

	stepStart = step0
	stepEnd = step1
	dStep = dstep

	filename = "optimalStep.dat"

	if sim:
		os.system(f"../compiled/optimalStep.exe {MCCs} {omega} {alphaStart} {alphaEnd} {dAlpha} {stepStart} {stepEnd} {dStep} {filename}")

	steps, alphas, percentAcceptAlphas = findPercentAccepts(filename)
	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(ncols=2, nrows=1,dpi = 125)
		ax = axes[0]
		ax2 = axes[1]
		ax.axhline(0.5, linestyle="dashed", c="k", zorder=10000, label=r"50% acceptance rate")
		ax.set_xlabel("Step size, $\delta$", fontsize=13)
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
		

		
		ax2.set_xlabel(rf"$\alpha$", fontsize=15)
		ax2.set_ylabel(r"$\delta'$", fontsize=15)
		optimalStep = findOptimalStep(steps, alphas, percentAcceptAlphas)
		alphas = alphas[1:]
		optimalStep = optimalStep[1:]
		for i,(alpha, oS) in enumerate(zip(alphas, optimalStep)):
			ax2.scatter(alpha, oS, color=colors[i])
		fitfunc = lambda x,m,b: m/x**(b)
		print(alphas)
		popt, pcov = curve_fit(fitfunc, alphas, optimalStep)
		print(popt)
		ax2.plot(alphas, fitfunc(alphas, *popt), color="k", alpha=0.7, lw=2,label=f"Fit: $\delta'(\\alpha) \\approx {round(popt[0],2)}"+r"\cdot \alpha^{-0.5}$")
		ax2.legend(loc='upper right', fancybox=True,shadow=True,fontsize=11)
		ax2.set_title(f"Step size, $\delta '$, yielding 50% acceptance rate. $\omega = {omega}$")
		plt.show()
