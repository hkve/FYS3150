def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os,sys
	from mpl_toolkits.mplot3d import Axes3D,art3d

	idx, MCCs, alpha, beta, E, EE, r12 = np.loadtxt(f"../data/betaGridMCCs8_100.dat", unpack=True)

	idx = np.argsort(idx)
	
	E = E[idx]
	alpha = alpha[idx]
	beta = beta[idx]
	alpha = np.unique(alpha)
	beta = np.unique(beta)
	print(beta)
	alpha, beta = np.meshgrid(alpha, beta)
	E = E.reshape(alpha.shape)
	minIdxVar = np.unravel_index(np.argmin(E, axis=None), E.shape, order="F")
	minIdxE = np.unravel_index(np.argmin(E, axis=None), E.shape, order="C")
	print(alpha[minIdxVar], beta[minIdxVar], E[minIdxE]) # This is alphaMin and betaMin

	
	fig = plt.figure()	
	ax = fig.gca(projection='3d')
	ax.scatter(alpha[minIdxVar], beta[minIdxVar],E[minIdxE], zorder=1000, s=20,label=rf"Minima: $\alpha = {round(alpha[minIdxVar],3)}, \beta = {round(beta[minIdxVar],3)}, \langle E \rangle = {round(E[minIdxE],4)}$", color="black")
	surf = ax.plot_surface(alpha, beta, E,cmap="jet",linewidth=8, antialiased=True,alpha=0.95)#"turbo"

	ax.set_xlabel(r"$\alpha$", fontsize=15)
	ax.set_ylabel(r"$\beta$", fontsize=15)
	ax.set_zlabel("$E$", fontsize=15)
	ax.set_title(r"Estimated ground state energy for varying $\alpha$ and $\beta$",fontsize=14)
	ax.legend(loc='upper center', fancybox=True,shadow=True,fontsize=12)
	plt.show()

main()
