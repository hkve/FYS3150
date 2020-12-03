def plotOptimalStep(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	from optimalStep import findPercentAccepts, findOptimalStep

	omega = 1
	alphaStart = 0.1
	alphaEnd = 2
	dAlpha = 0.1

	stepStart = 0.1
	stepEnd = 5
	dStep = 0.3

	filename = "optimalStep1.dat"

	if sim:
		os.system(f"../compiled/optimalStep.exe {6} {omega} {alphaStart} {alphaEnd} {dAlpha} {stepStart} {stepEnd} {dStep} {filename}")

	steps, alphas, percentAcceptAlphas = findPercentAccepts(filename)

	fig, ax = plt.subplots()
	ax.axhline(0.5, linestyle="dashed", c="k")
	ax.set_xlabel("Step length", fontsize=13)
	ax.set_ylabel("Accepted", fontsize=13)
	
	for alpha, percentAccepts in zip(alphas, percentAcceptAlphas):		
		ax.plot(steps, percentAccepts, label=rf"$\alpha=${alpha:.2f}")
	

	ax.legend()
	plt.show()

	fig, ax = plt.subplots()
	ax.set_xlabel(rf"$\alpha$ ($\omega=${omega:.3f})", fontsize=13)
	ax.set_ylabel("Step closest to 50% accepts", fontsize=13)
	optimalStep = findOptimalStep(steps, alphas, percentAcceptAlphas)
	
	ax.scatter(alphas, optimalStep)
	plt.show()

plotOptimalStep(sim = True)