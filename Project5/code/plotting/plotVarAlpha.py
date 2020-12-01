def varAlpha(sim = False): 
	import numpy as np
	import matplotlib.pyplot as plt
	import os

	MCCs = 7
	alphaStart = 0.5
	alphaEnd = 1.5
	dAlpha = 0.1

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1

	alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)

	if sim:
		os.chdir("../VMC/compiled")
		os.system(f"./variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} noninteractive")
		os.chdir("../../plotting")

	fig, axes = plt.subplots(nrows=1, ncols=2)
	axes[0].set(xscale="log")
	axes[1].set(xscale="log")
	for i, alpha in enumerate(alphas):
		data = np.loadtxt(f"../data/varAlphaNoninteractive{i}.dat")
		cyles, E, E2 = data[:,0], data[:,1], data[:,2]

		varE = E2-E*E
		skip = 0
		axes[0].plot(cyles[skip:], E[skip:], label=fr"$\alpha$ = {alpha:.1f}, E = {E[-1]:.2f}")
		axes[1].plot(cyles[skip:], varE[skip:], label=fr"$\alpha$ = {alpha:.1f} $\sigma_E^2$ = {varE[-1]:.2f}")

	axes[0].legend()
	axes[1].legend()
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
		os.chdir("../VMC/compiled")
		os.system(f"./variateAlpha.exe {7} {alphaStart} {alphaEnd} {dAlpha} interactive")
		os.chdir("../../plotting")

	fig, axes = plt.subplots(nrows=1, ncols=2)
	axes[0].set(xscale="log")
	axes[1].set(xscale="log")
	for i, alpha in enumerate(alphas):
		data = np.loadtxt(f"../data/varAlphaInteractive{i}.dat")
		cyles, E, E2 = data[:,0], data[:,1], data[:,2]

		varE = E2-E*E
		skip = 0
		axes[0].plot(cyles[skip:], E[skip:], label=fr"$\alpha$ = {alpha:.2f}, E = {E[-1]:.2f}")
		axes[1].plot(cyles[skip:], varE[skip:], label=fr"$\alpha$ = {alpha:.2f} $\sigma_E^2$ = {varE[-1]:.2f}")

	axes[0].legend()
	axes[1].legend()
	plt.show()
