def findPercentAccepts(filename):
	# I eat datafile nomnom and make arrays
	import numpy as np

	data = np.loadtxt(f"../data/{filename}")

	MCCs, alpha, step, accepts = data[:,0], data[:,1], data[:,2], data[:,3]

	alphas = np.unique(alpha)
	steps = np.unique(step)
	
	percentAccepts = accepts/MCCs
	percentAcceptsAlpha = []
	for a in alphas:
		dataIdx = np.in1d(alpha, round(a,10))
		
		percentAcceptsAlpha.append(percentAccepts[dataIdx])
	
	return steps, alphas, np.array(percentAcceptsAlpha)

def findOptimalStep(steps, alphas, percentAccepts):
	# ahh *smatt smatt* good datafile, give me arrays and i make more array
	import numpy as np

	bestStep = []
	for alpha, percentAccept in zip(alphas, percentAccepts):
		bestStepIdx = np.argmin(abs(percentAccept-0.5))
		bestStep.append(steps[bestStepIdx])

	return np.array(bestStep)
