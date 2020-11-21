def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run

	MCCs = int(1e6)
	T_ = [1,2.4]
	initSpin_ = [0,1]

	if sim: 
		for T in T_:
			for initSpin in initSpin_:
				filename = f"20x20_{T}_{initSpin}.dat"	

				run(f"../cpp/20x20_lattice.out {T} {MCCs} {initSpin} {filename}".split())

	filenames = []
	for T in T_:
		for initSpin in initSpin_:
			filenames.append(f"20x20_{T}_{initSpin}.dat")

	fig, ax = plt.subplots(nrows=1, ncols=2)
	for filename in filenames:
		data = np.loadtxt(f"../data/{filename}")
		E, Mabs, MCCs = data[:,0], data[:,1], data[:,2]

		ax[0].plot(MCCs, E/20**2)
		ax[1].plot(MCCs, Mabs/20**2)

	plt.show()

main(sim = True)