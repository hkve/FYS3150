def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os,sys
	from mpl_toolkits.mplot3d import Axes3D,art3d
	from matplotlib.patches import Circle

	MCCs = 6
	omega = 1

	# These values gave pretty plots, feel free to test different
	alphaMin = 0.8 
	alphaMax = 1.4
	N_alpha = 250
	betaMin = 0
	betaMax = 1
	N_beta = 250
	filename = f"variateBeta_grid_{betaMin}-{betaMax}-{N_beta}_{alphaMin}-{alphaMax}-{N_alpha}_{MCCs}-{omega}.data"

	if sim:
		if filename in os.listdir("../data/"):
			os.system(f"rm ../data/{filename}")
		print(f"../compiled/variateBeta_grid.exe {MCCs} {omega} {alphaMin} {alphaMax} {N_alpha} {betaMin} {betaMax} {N_beta} {filename}")
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
	ax.scatter(0.956, 0.333,3.730, zorder=1000, s=20,label=f"Minima: $\\alpha = 0.956, \\beta = 0.333, E = 3.730$", color="black")
	
	surf = ax.plot_surface(alpha, beta, E,cmap="turbo",linewidth=8, antialiased=True,alpha=0.95)#"turbo"
	
	ax.set_xlabel(r"$\alpha$", fontsize=15)
	ax.set_ylabel(r"$\beta$", fontsize=15)
	ax.set_zlabel("$E$", fontsize=15)
	ax.set_title(r"Estimated ground state energy for varying $\alpha$ and $\beta$",fontsize=14)
	ax.legend(loc='upper center', fancybox=True,shadow=True,fontsize=12)
	plt.show()
	
