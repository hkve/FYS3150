def varBeta(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from mpl_toolkits.mplot3d import Axes3D

	MCCs = 6
	omega = 1

	# These values gave pretty plots, feel free to test different
	alphaMin = 0.8 
	alphaMax = 1.4
	N_alpha = 20
	betaMin = 0
	betaMax = 1
	N_beta = 20
	filename = "variateBeta_grid.dat"

	if sim:
		if filename in os.listdir("../data/"):
			os.system(f"rm ../data/{filename}")
		os.system(f"../compiled/variateBeta_grid.exe {MCCs} {omega} {alphaMin} {alphaMax} {N_alpha} {betaMin} {betaMax} {N_beta} {filename}")

	MCCs, alpha, beta, E, EE, varE = np.loadtxt(f"../data/{filename}", unpack=True)

	alpha = np.unique(alpha)
	beta = np.unique(beta)
	
	alpha, beta = np.meshgrid(alpha, beta)
	
	E = E.reshape(alpha.shape)
	
	minIdx = np.unravel_index(np.argmin(E, axis=None), E.shape)
	print(alpha[minIdx], beta[minIdx]) # This is alphaMin and betaMin

	fig = plt.figure()	
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(alpha, beta, E,cmap="gnuplot")
	ax.set_xlabel(r"$\alpha$", fontsize=13)
	ax.set_ylabel(r"$\beta$", fontsize=13)
	ax.set_zlabel("$E$", fontsize=13)
	plt.show()
	
varBeta(sim = False)