def varAlpha(sim = False): 
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	#MCCs = 7
	alphaStart = 0.5
	alphaEnd = 1.5
	dAlpha = 0.01

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1

	fig, axes = plt.subplots(nrows=1, ncols=2)
	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)
	E_ = []
	r12_ = []
	for mode in ["interactive", "noninteractive"]:
		if sim:
			os.system(f"../compiled/variateAlpha.exe {5} {alphaStart} {alphaEnd} {dAlpha} {mode} NoSave")

	#axes.set(xscale="log")
	#axes[1].set(xscale="log")
	
		cycles, beta, E, E2, r12 = np.loadtxt(f"../data/alpha{mode}.dat", unpack=True, usecols=(0,2,3,4,5))
		#	file << MCCs << " " << alpha << " " << beta << " " << E << " " << EE << " " << r12 <<endl;
		E_.append(E)
		r12_.append(r12)
		#cyles, beta E, E2 = data[:,0], data[:,1], data[:,2]
		#print(cycles, beta, E, E2, r12)
		
		
		axes[0].plot(alphas, E,label=mode) #,label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
		#axes[1].plot(cycles[i], 	varE)#, label=fr"$\alpha$ = {alpha:.1f} $\sigma_E^2$ = {varE[-1]:.2f}")
		axes[1].plot(alphas, r12, label=mode)
		#axes[1].plot(alphas, 1/r12, label=fr"$1/r_{12}000$")
		#axes[1].legend()
	#axes[1].plot(alphas, E_[0]-E_[1], label="E2-E1")

	#axes[1].plot(alphas, 1/r12, label=fr"$1/r_{12}1111$")
	#EL2 = lambda b:
	#b=0
	#EL2 = E_[0] + 1/(2*(1+b*r12)**2) *(alphas*r12 -1/(2*(1+b*r12)**2) -2/r12+2*b/(1+b*r12))
	#axes[1].plot(alphas, E_[0]- E_[1], label="E0-E1")

	#axes[1].plot(alphas, 1/r12,label=fr"1/r12")

	#print(np.mean(EL2- E_[1] - 1/r12), np.mean(EL2))
	axes[1].legend()
	axes[0].legend()
	plt.show()

def varAlphaInteractive(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	MCCs = 7
	alphaStart = 0.85
	alphaEnd = 1.15
	dAlpha = 0.05

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1

	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)

	if sim:
		os.system(f"../compiled/variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} interactive")

	fig, axes = plt.subplots(nrows=1, ncols=2)
	axes[0].set(xscale="log")
	axes[1].set(xscale="log")
	
	data = np.loadtxt(f"../data/varAlphaInteractive.dat")
	cyles, E, E2 = data[:,0], data[:,1], data[:,2]

	varE = E2-E*E
	skip = 0
	#axes[0].plot(cyles[skip:], E[skip:], label=fr"$\alpha$ = {alpha:.2f}, E = {E[-1]:.2f}")
	#axes[1].plot(cyles[skip:], varE[skip:], label=fr"$\alpha$ = {alpha:.2f} $\sigma_E^2$ = {varE[-1]:.2f}")

	axes[0].legend()
	axes[1].legend()
	plt.show()

varAlpha(sim = False)