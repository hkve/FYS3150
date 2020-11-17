def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import sys

	L = 20
	orientation = 0 # random
	T_ = [1, 2, 3]
	N = 10
	MCCs = np.logspace(1, 4, N, dtype=int, endpoint=True)
	filenames = [f"acceptedFlips_{T}.out" for T in T_]
	for T, filename in zip(T_, filenames):
		for MCC in MCCs:
			if sim:
				print(T, MCC)
				run(f"../cpp/acceptedFlips.out {T} {L} {MCC} {orientation} {filename}".split())
			else:
				break

	 
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots()
		ax.set(xscale="log")
		for T, filename in zip(T_, filenames):
			data = np.loadtxt(f"../data/{filename}")
	
			acceptedFlips, attempedFlips, ratio, L, MCCs = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4] 
			ax.scatter(MCCs, acceptedFlips, label=f"T={T}")
		plt.legend()
		plt.show()
		
