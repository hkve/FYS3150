def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import sys
	from colour import Color
	import matplotlib as mpl

	L = 20
	orientation = 0 # random
	T_ = [1, 2, 3]
	T_ = np.linspace(1,2.4, 5, endpoint=True)

	N = 100
	MCCs = np.logspace(1, 5, N, dtype=int, endpoint=True)
	filenames = [f"acceptedFlips_{T}.out" for T in T_]
	for T, filename in zip(T_, filenames):
		for MCC in MCCs:
			if sim:
				sys.exit()
				print(round(T,3), MCC/MCCs[-1])
				run(f"../cpp/acceptedFlips.out {T} {L} {MCC} {orientation} {filename}".split())
			else:
				break

	
	with sns.axes_style("darkgrid"):
		fig, ax = plt.subplots(nrows=1, ncols=2,dpi=150)
		ax[0].set(xscale="log", yscale="log")
		ax[1].set(xscale="log")

		colors = ["red", "orange", "green", "cyan", "blue"]
		R = np.zeros((len(T_),len(MCCs[-20:])))
		#colors2 = [Color("red", luminance=l), Color("orange", luminance=l),Color("green", luminance=l),Color("blue", luminance=l),Color("purple", luminance=l)]
		#print(colors2)
		#colors = list(Color("red", luminance=0.4).range_to(Color("blue"),len(T_))) 
		#colors = [color.get_rgb() for color in colors]
		for i,(T, filename) in enumerate(zip(T_[::-1], filenames[::-1])):
			data = np.loadtxt(f"../data/{filename}")
			print(MCCs[-20])
			acceptedFlips, attempedFlips, ratio, L, MCCs = data[:,0],data[:,1],data[:,2],data[:,3],data[:,4] 
			print(ratio)
			m,b = np.polyfit(np.log10(MCCs[-20:]), np.log10(acceptedFlips[-20:]), 1)
			m2,b2 = np.polyfit(MCCs[-20:], acceptedFlips[-20:], 1)

			
			ax[0].plot(MCCs, acceptedFlips, alpha=0.7, lw=3, label="T=%.2f" %(T), color=colors[i])
			#R[i] = acceptedFlips[-20:]/(L**2*MCCs[-20:])
			ax[1].plot(MCCs,ratio, color=colors[i],alpha=0.7, lw=3)
		# print(R.shape)
		# for a in range(R.shape[0]):
		# 	#print(R[a])
		# 	ax.plot(T_, R[:,a], label="a")
		# cmap = mpl.colors.ListedColormap(colors)
		# norm = mpl.colors.Normalize(vmin=T_[0],vmax=T_[-1])
		# cbar = ax[1].figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[1], orientation='vertical',ticks=T_)
		ax[0].set_xlabel("MCCs")
		ax[0].set_ylabel("Accepted flips")
		ax[0].set_title("Number of accepted flips")
		ax[1].set_title("Acceptance rate, $R$")
		ax[1].set_xlabel("MCCs")
		ax[1].set_ylabel("$R$")
		lines, labels = fig.axes[0].get_legend_handles_labels()
		fig.legend(lines,labels,loc='lower center', ncol=5, fancybox=True, shadow=True,fontsize=10)
		plt.show()
		
