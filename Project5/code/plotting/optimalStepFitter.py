import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import seaborn as sns
w = [0.01, 0.1, 0.25, 0.5, 0.75, 1]
m = [13.8, 4.387, 2.77, 1.96, 1.60, 1.39]
fitfunc = lambda x,m: m/x**(0.5)
w_hires = np.linspace(w[0],w[-1], 10000)
popt, pcov = curve_fit(fitfunc, w, m)
print(popt)
with sns.axes_style("darkgrid"):
	fig, ax = plt.subplots(ncols=1, nrows=1,dpi = 140)
	ax.plot(w,m, "-o", label="Dataset")
	ax.plot(w_hires, fitfunc(w_hires, *popt), color="k", alpha=0.7, lw=2, label=r"Fit: $m(\omega) = 1.381\cdot \omega^{-0.5}$")
	ax.legend(fancybox=True, shadow=True)
	ax.set_xlabel(r"$\omega$", fontsize=14)
	ax.set_ylabel(r"$m$", fontsize=14)
	ax.set_title(r"Best fit of $m$ in $\delta'(\alpha,\omega) = m(\omega)\alpha^{-0.5}$", fontsize=11)
	plt.show()