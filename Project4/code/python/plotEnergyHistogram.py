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

	
	
	colors = ["b", "r"]

	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(ncols=2, nrows=1, dpi=120)

		for ax_, T, filename, color in zip(ax, T_, filenames, colors):
			data = np.loadtxt(f"../data/{filename}")
			
			
			weights = np.ones_like(data)/len(data)
			var = "%.2e" %np.var(data)
			mean = str(round(np.mean(data),2))
			Elabel = f"$\\langle E \\rangle = {mean}$"
			varlabel = f"$\\sigma^2 = {var}$"
			binwidth = 0.01
			ax_.hist(data, bins=np.arange(min(data), max(data) + binwidth, binwidth), weights=weights, label=varlabel, facecolor=color, edgecolor=color)
			ax_.plot((),(), color=color, label=Elabel)
			ax_.set_title(f"$T = {T}$", fontsize=13)
			leg = ax_.legend(handlelength=0, handletextpad=0, loc='upper right', fancybox=True, shadow=True,fontsize=12)
			for item in leg.legendHandles:
				item.set_visible(False)
			print(f"T = {T} E_mean = {np.mean(data)}, sigma^2 = {np.var(data)}")
			fontsize=13
			ax_.set_xlabel("$E$", fontsize=fontsize)
			ax_.set_ylabel("$P(E)$", fontsize=fontsize)
		plt.show()
