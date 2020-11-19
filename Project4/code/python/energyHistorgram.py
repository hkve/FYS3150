def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	import seaborn as sns
	from subprocess import run
	import os

	T_ = [1, 2.4]
	MCCs = int(1e6)
	L = 20
	stableAfter = int(10**4.2)
	filenames = [f"20x20_{T}_energy_count.dat" for T in T_]

	for T, filename in zip(T_, filenames):
		if sim:
			if filename in os.listdir("../data/"):
				os.system(f"rm ../data/{filename}")
			run(f"../cpp/20x20_energy_count.out {T} {stableAfter} {MCCs} {filename}".split())

	bins = [5, 59] # Disse bin-sizene ga ish fine plots, bare å eksperimentere
	# Valgte log siden jeg ville plotte begge i samme, kan evt å ha 2 y-akser 
	fig, ax = plt.subplots()
	ax.set(yscale="log")
	for T, filename, Bin in zip(T_, filenames, bins):
		data = np.loadtxt(f"../data/{filename}")
		weights = np.ones_like(data)/len(data)
		ax.hist(data, bins=Bin, weights=weights, alpha=0.75, label=f"{T}")
		print(f"T = {T} E_mean = {np.mean(data)}, sigma = {np.std(data)}")
		ax.set(xlabel="E", ylabel="P(E)")
	ax.legend()
	plt.show()
	
main(sim = False)