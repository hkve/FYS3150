def main(sim = False, comp = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import os

	T_ = [1]
	N = 200
	MCCs = np.logspace(1, 6, N, dtype=int, endpoint=True)
	
	filenames = [f"2x2_{T}.out" for T in T_]
	if comp: 
		os.chdir("../cpp/")
		os.system("make 2x2")
		os.chdir("../python/")
	for T, filename in zip(T_, filenames):
		if sim:
			os.system(f"rm ../data/{filename}")
		for MCC in MCCs:
			if sim:
				print(T, MCC)
				run(f"../cpp/2x2_lattice.out {T} {T} {0} {MCC} {filename}".split())
			else:
				break


	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=2, ncols=2, dpi=100)
		axes[0,0].set(xscale="log")
		axes[0,1].set(xscale="log")
		axes[1,0].set(xscale="log")
		axes[1,1].set(xscale="log")
		for filename in filenames:
			data = np.loadtxt(f"../data/{filename}")
			E, M, E2, M2, Mabs, varE, varM, T, MCCs = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
			L = 2

			E = E/L**2
			Mabs = Mabs/L**2
			Z = 12 + 4*np.cosh(8/T)
			CV = varE/T**2/L**2
			X = varM/T/L**2
			Cv_ana = 1024/(T*Z)**2 *(3*np.cosh(8/T)+1)/L**2
			E_ana = -32/Z * np.sinh(8/T) /L**2
			Mabs_ana = 8/Z * (2+np.exp(8/T))
			X_ana = (32/(T*Z) * (1+np.exp(8/T)) - Mabs_ana**2)**2/(T*L**2)
			Mabs_ana /= L**2

			markersize = 50
			markeralpha = 0.6
			anastyle = "solid"
			anasize = 3
			anaalpha = 0.6
			axes[0,0].scatter((),(), label="Data", marker="+", c="k",s=80)

			axes[0,0].scatter(MCCs, E,  marker="+", c="k",s=markersize, alpha=markeralpha)
			axes[0,0].plot([MCCs[0],MCCs[-1]], [E_ana[0], E_ana[-1]], c="r", linestyle=anastyle, lw = anasize, label=f"Analytical", alpha=anaalpha)
			axes[0,1].scatter(MCCs, Mabs, marker="+", c="k",s=markersize, alpha=markeralpha)
			axes[0,1].plot([MCCs[0],MCCs[-1]], [Mabs_ana[0], Mabs_ana[-1]], c="r", linestyle=anastyle, lw = anasize, alpha=anaalpha)
			axes[1,0].scatter(MCCs, CV, marker="+", c="k",s=markersize, alpha=markeralpha)
			axes[1,0].plot([MCCs[0], MCCs[-1]], [Cv_ana[0],Cv_ana[-1]], c="r",linestyle=anastyle, lw = anasize, alpha=anaalpha)
			axes[1,1].scatter(MCCs, X, marker="+", c="k", s=markersize, alpha=markeralpha)
			axes[1,1].plot([MCCs[0], MCCs[-1]], [X_ana[0],X_ana[-1]], c="r",linestyle=anastyle, lw = anasize, alpha=anaalpha)
			fontsize=14
			axes[0,0].set_ylabel(r"$\langle E \rangle$", fontsize=fontsize)
			axes[0,1].set_ylabel(r"$\langle | M| \rangle$", fontsize=fontsize)
			axes[1,0].set_ylabel("$C_v$", fontsize=fontsize)
			axes[1,1].set_ylabel("$\chi$", fontsize=fontsize)
			axes[1,0].set_xlabel("MCCs", fontsize=fontsize)
			axes[1,1].set_xlabel("MCCs", fontsize=fontsize)
		
		axes[0,0].set_xticks([10**1, 10**2, 10**3, 10**4,10**5,10**6])
		axes[0,1].set_xticks([10**1, 10**2, 10**3, 10**4,10**5,10**6])
		axes[1,0].set_xticks([10**1, 10**2, 10**3, 10**4,10**5,10**6])
		axes[1,1].set_xticks([10**1, 10**2, 10**3, 10**4,10**5,10**6])

		axes[0,1].set_ylim([0.85,1.01])
		axes[0,1].set_yticks([0.85, 0.9, 0.95, 1.0])
		lines, labels = axes[0,0].get_legend_handles_labels()
		fig.legend(lines,labels,loc='upper center', ncol=2, fancybox=True, shadow=True,fontsize=14)
		plt.show()
	