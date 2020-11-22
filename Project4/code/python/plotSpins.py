import numpy as np
import os 
import matplotlib.pyplot as plt 

def main(sim = False):
	L = 500
	T_ = [0.1, 500]
	MCCs = int(200)
	stableMCCs = int(1)
	filenames = [f"spins_L{L}_T{T}_MCCs{MCCs}_stableMCCS{stableMCCs}.dat" for T in T_]
	datas = []
	for T,filename in zip(T_,filenames):
		if sim:
			if os.path.isfile(f"../data/{filename}"):
				os.system(f"rm ../data/{filename}")
			os.system(f"../cpp/plot_spins.out {L} {T} {MCCs} {stableMCCs} {filename}")
		datas.append(np.loadtxt(f"../data/{filename}"))

	

	fig, ax = plt.subplots(ncols=2, nrows=1, dpi=220)

	titles = [f"$T = {T}$" for T in T_]
	for i,(ax_,title) in enumerate(zip(ax,titles)):
		ax_.imshow(datas[i].reshape(L,L), interpolation="nearest",cmap='Greys')
		ax_.set_yticks([])
		ax_.set_xticks([])
		ax_.set_title(title)
	
	plt.show()
main(sim=False)
