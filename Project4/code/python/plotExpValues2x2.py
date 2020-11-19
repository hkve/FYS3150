def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import seaborn as sns
	from subprocess import run
	
	Tstart = 1
	Tend = 4
	dt = 0.1
	MCCs = int(1e6)
	filename = "2x2EXP.out"

	if sim:
		#print("--sim called. however, no default function call for plotExpValues2x2.py has been defined")
		run(f"../cpp/2x2_lattice.out {Tstart} {Tend} {dt} {MCCs} {filename}".split())
			
	data = np.loadtxt("../data/2x2EXP.out")
	L = 2
	E, M, E2, M2, Mabs, varE, varM, T = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]
	CV = varE/T**2
	X = varM/T
	Z = 12 + 4*np.cosh(8/T)
	E_ana = -32/Z * np.sinh(8/T) 
	Mabs_ana = 8/Z * (2+np.exp(8/T))
	Cv_ana = 1024/(T*Z)**2 *(3*np.cosh(8/T)+1)
	#X_ana = 32/(T*Z) * (1+np.exp(8/T))
	X_ana = (32/Z * (1+np.exp(8/T)) - Mabs_ana**2)/T
	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=2,ncols=2, dpi=100)
		markersize =120
		markeralpha = 0.9
		anastyle = "solid"
		anasize = 4
		anaalpha = 1
		axes[0,0].scatter(T, E/L**2,label="Data", marker=".", s=markersize, c="k", alpha=markeralpha,zorder=100)
		axes[0,1].scatter(T, Mabs/L**2, marker=".", s=markersize, c="k", alpha=markeralpha,zorder=100)
		axes[1,0].scatter(T, CV/L**2, marker=".",s=markersize, c="k", alpha=markeralpha,zorder=100)
		axes[1,1].scatter(T, X/L**2,marker=".", s=markersize, c="k", alpha=markeralpha,zorder=100)

		axes[0,0].plot(T, E_ana/L**2,label="Analytical", c="r", lw=anasize, alpha=anaalpha)
		axes[0,1].plot(T, Mabs_ana/L**2, c="r", lw=anasize, alpha=anaalpha)
		axes[1,0].plot(T, Cv_ana/L**2, c="r", lw=anasize, alpha=anaalpha)
		axes[1,1].plot(T, X_ana/L**2, c="r", lw=anasize, alpha=anaalpha)

		fontsize=14
		axes[0,0].set_ylabel(r"$\langle E \rangle$", fontsize=fontsize)
		axes[0,1].set_ylabel(r"$\langle |M| \rangle$", fontsize=fontsize)
		axes[1,0].set_ylabel("$C_v$", fontsize=fontsize)
		axes[1,1].set_ylabel("$\chi$", fontsize=fontsize)
		axes[1,0].set_xlabel("$T$  $[ kT$ / $J]$", fontsize=fontsize)
		axes[1,1].set_xlabel("$T$  $[ kT$ / $J]$", fontsize=13)
	lines, labels = axes[0,0].get_legend_handles_labels()
	fig.legend(lines,labels,loc='upper center', ncol=2, fancybox=True, shadow=True,fontsize=14)
	plt.show()

main(sim=False)