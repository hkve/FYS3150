def main(sim = False, comp = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import os

	L = 20
	MCCs = int(1e6) # FHODSFHIODSFHDSFHOHFODHFOHDOJFHOSD
	T_ = [1,2.4]
	initSpin_ = [0,1]
	filenames = []
	
	if comp:
		os.chdir("../cpp/")
		os.system("make 20x20")
		os.chdir("../python/")
	for T in T_:
		for initSpin in initSpin_:
			filename = f"20x20_{T}_{initSpin}.dat"
			filenames.append(filename)
			if sim:
				try:
					os.system(f"rm ../data/{filename}")
				finally:		
					run(f"../cpp/20x20_lattice.out {T} {MCCs} {initSpin} {filename}".split())
					print(f"Done T = {T}, initSpin = {initSpin}")
	
	spinLabel = ["Random", "Uniform"]
	Tcolors = ["blue", "red"]
	spinMarker = ["solid", (0,(0.1,1.6))]
	spinAlpha=[0.6, 0.8]
	with sns.axes_style("darkgrid"):
		fig = plt.figure(dpi=120)
		axes = fig.subplots(nrows=1,ncols=2, sharex=True)
		for i,initSpin,filename in zip([0,0,1,1],[0,1,0,1],filenames):
			data = np.loadtxt(f"../data/{filename}")
			E, Mabs, MCCs = data[:,0], data[:,1], data[:,2]

			#axes[0].plot(MCCs, E/20**2, color=)
			#axes[1].plot(MCCs, Mabs/20**2)
			axes[0].plot(MCCs, E/L**2, linestyle=spinMarker[initSpin] , dash_capstyle = 'round',lw=3.5, color = Tcolors[i], alpha=spinAlpha[initSpin])#,label=f"$T = {T}$, {spinLabel[initSpin]} spins")
			axes[0].plot(0, linestyle=spinMarker[initSpin] ,lw=3.5, dash_capstyle = 'round',color = Tcolors[i],label=f"$T = {T_[i]}$, {spinLabel[initSpin]} spins")
			axes[1].plot(MCCs, Mabs/L**2, linestyle=spinMarker[initSpin], lw=3.5, color=Tcolors[i], dash_capstyle = 'round', alpha=spinAlpha[initSpin])#, label=f"$T = {T}$, {spinLabel[initSpin]} spins")
		#         foo += 1
		# for i,T in enumerate(T_):
		#     for j,initSpin in enumerate(initSpin_):
		#         axes[0].plot(MCCs, E_[foo]/L**2, linestyle=spinMarker[initSpin] , dash_capstyle = 'round',lw=3.5, color = Tcolors[i], alpha=spinAlpha[initSpin])#,label=f"$T = {T}$, {spinLabel[initSpin]} spins")
		#         axes[0].plot(0, linestyle=spinMarker[initSpin] ,lw=3.5, dash_capstyle = 'round',color = Tcolors[i],label=f"$T = {T}$, {spinLabel[initSpin]} spins")

		 #         axes[1].plot(MCCs, Mabs_[foo]/L**2, linestyle=spinMarker[initSpin], lw=3.5, color=Tcolors[i], dash_capstyle = 'round', alpha=spinAlpha[initSpin])#, label=f"$T = {T}$, {spinLabel[initSpin]} spins")
		#         foo += 1
		axes[0].set_ylim(-2.1, -0.5)
		axes[1].set_ylim(-0.1,1.1)
		lines, labels = fig.axes[0].get_legend_handles_labels()
		fig.legend(lines,labels,loc='upper center', ncol=2, fancybox=True, shadow=True,fontsize=10)#, bbox_to_anchor=(0.5, 1.16),ncol=2, fancybox=True, shadow=True, fontsize=9)
		axes[1].set_xlabel("MCCs", fontsize=12)
		axes[0].set_xlabel("MCCs", fontsize=12)

		axes[0].set_ylabel("$\\langle E \\rangle / L^2$", fontsize=13)
		axes[1].set_ylabel("$ \\left <|M|\\right > / L^2 $", fontsize=13)
	# print("Plot not in log(x)-y format, press 'k' in plot window to change format")
	# plt.show()
	# fig, ax = plt.subplots(nrows=1, ncols=2)
	# for filename in filenames:
	# 	data = np.loadtxt(f"../data/{filename}")
	# 	E, Mabs, MCCs = data[:,0], data[:,1], data[:,2]

	# 	ax[0].plot(MCCs, E/20**2)
	# 	ax[1].plot(MCCs, Mabs/20**2)

	plt.show()
