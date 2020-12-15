def varBetaParallel(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os,sys
	from mpl_toolkits.mplot3d import Axes3D,art3d

	idx, MCCs, alpha, beta, E, EE, r12 = np.loadtxt(f"../data/dump.dat", unpack=True)

	idx = np.argsort(idx)
	
	E = E[idx]
	alpha = alpha[idx]
	beta = beta[idx]
	alpha = np.unique(alpha)
	beta = np.unique(beta)
	
	alpha, beta = np.meshgrid(alpha, beta)
	
	E = E.reshape(alpha.shape)
	minIdx = np.unravel_index(np.argmin(E, axis=None), E.shape)
	print(alpha[minIdx], beta[minIdx]) # This is alphaMin and betaMin
	
	fig = plt.figure()	
	ax = fig.gca(projection='3d')
	ax.scatter(alpha[minIdx], beta[minIdx],E[minIdx], zorder=1000, s=20,label=f"Minima: $\\alpha = {round(alpha[minIdx],2)}, \\beta = {round(beta[minIdx],2)}, E = {round(E[minIdx],2)}$", color="black")
	surf = ax.plot_surface(alpha, beta, E,cmap="jet",linewidth=8, antialiased=True,alpha=0.95)#"turbo"

	ax.set_xlabel(r"$\alpha$", fontsize=15)
	ax.set_ylabel(r"$\beta$", fontsize=15)
	ax.set_zlabel("$E$", fontsize=15)
	ax.set_title(r"Estimated ground state energy for varying $\alpha$ and $\beta$",fontsize=14)
	ax.legend(loc='upper center', fancybox=True,shadow=True,fontsize=12)
	plt.show()

varBetaParallel()