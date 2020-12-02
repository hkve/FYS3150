def plotOptimalStep(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt

	alphaStart = 0.7
	alphaEnd = 1.3
	dAlpha = 0.1

	stepStart = 0.5
	stepEnd = 20
	dStep = 0.5

	N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
	N_alpha = int(N_alpha) + 1
	N_step = round(((stepEnd-stepStart)/dStep),0)
	N_step = int(N_step) + 1
	
	alphas = np.linspace(alphaStart, alphaEnd, N_alpha)
	steps = np.linspace(stepStart, stepEnd, N_step)

	data = np.loadtxt("../data/optimalStep1.dat")

	MCCs, alpha, step, accepts = data[:,0], data[:,1], data[:,2], data[:,3]

	percentAccepts = accepts/MCCs

	fig, ax = plt.subplots()
	
	ax.axhline(0.5, linestyle="dashed", c="k")
	ax.set_xlabel("Step length", fontsize=13)
	ax.set_ylabel("Accepted", fontsize=13)
	for a in alphas:
		dataIdx = np.in1d(alpha, round(a,10))
		ax.plot(steps, percentAccepts[dataIdx], label=rf"$\alpha=${a:.2f}")

	ax.legend()
	plt.show()

	
plotOptimalStep()