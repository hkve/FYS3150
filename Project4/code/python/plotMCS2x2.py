def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import sys

	MCCs_start = 10
	MCCs_end = 10000
	T_ = [1,3]
	N = 200
	MCCs = np.logspace(1, 6, N, dtype=int)
	
	filenames = [f"2x2_{T}.out" for T in T_]
	for T, filename in zip(T_, filenames):
		for MCC in MCCs:
			if sim:
				print(T, MCC)
				run(f"../cpp/2x2_lattice.out {T} {T} {0} {MCC} {filename}".split())
			else:
				break


	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=1, ncols=2)
		axes[0].set(xscale="log")
		axes[1].set(xscale="log")
		for filename in filenames:
			data = np.loadtxt("../data/" + filename)

			T = data[0,-2]
			MCCs = data[:,-1]
			E, Mabs = data[:,0], data[:,4]
			
			Z = 12 + 4*np.cosh(8/T)
			E_ana = -32/Z * np.sinh(8/T) 
			Mabs_ana = 8/Z * (2+np.exp(8/T))

			E_rel_error = np.abs((E-E_ana)/E_ana)
			Mabs_rel_error = np.abs((Mabs-Mabs_ana)/Mabs_ana)
			
			#E = E/np.linalg.norm(E)
			#Mabs = Mabs/np.max(Mabs)
			
			axes[0].scatter(MCCs, E, label=f"{T}")
			axes[0].plot([MCCs[0],MCCs[-1]], [E_ana, E_ana], c="k", linestyle="dashed")
			axes[1].scatter(MCCs, Mabs, label=f"{T}")
			axes[1].plot([MCCs[0],MCCs[-1]], [Mabs_ana, Mabs_ana], c="k", linestyle="dashed")
		
		lines, labels = fig.axes[0].get_legend_handles_labels()
		fig.legend(lines,labels,loc='upper center', ncol=2, fancybox=True, shadow=True,fontsize=13)
		plt.show()
	