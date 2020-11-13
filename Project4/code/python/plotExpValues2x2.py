def main():
	import numpy as np
	import matplotlib.pyplot as plt
	import seaborn as sns

	data = np.loadtxt("../data/expValues.out")

	E, M, E2, M2, Mabs, varE, varM, T = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]
	CV = varE/T**2
	X = varM/T
	Z = 12 + 4*np.cosh(8/T)
	E_ana = -32/Z * np.sinh(8/T) 
	Mabs_ana = 8/Z * (2+np.exp(8/T))
	Cv_ana = 1024/(T*Z)**2 *(3*np.cosh(8/T)+1)
	X_ana = 32/(T*Z) * (1+np.exp(8/T))

	with sns.axes_style("darkgrid"):
		fig, axes = plt.subplots(nrows=2,ncols=2)
		
		axes[0,0].scatter(T, E)
		axes[0,1].scatter(T, Mabs)
		axes[1,0].scatter(T, CV)
		axes[1,1].scatter(T, X)

		axes[0,0].plot(T, E_ana, c="r")
		axes[0,1].plot(T, Mabs_ana, c="r")
		axes[1,0].plot(T, Cv_ana, c="r")
		axes[1,1].plot(T, X_ana, c="r")

		fontsize=14
		axes[0,0].set_ylabel(r"$\langle E \rangle$", fontsize=fontsize)
		axes[0,1].set_ylabel(r"$\langle M \rangle$", fontsize=fontsize)
		axes[1,0].set_ylabel("$C_v$", fontsize=fontsize)
		axes[1,1].set_ylabel("$\chi$", fontsize=fontsize)
		axes[1,0].set_xlabel("T", fontsize=fontsize)
		axes[1,1].set_xlabel("T", fontsize=13)
	plt.show()