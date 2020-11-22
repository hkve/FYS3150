
def get_handle_lists(l):
	tree = l._legend_box.get_children()[1]

	for column in tree.get_children():
		for row in column.get_children():
			yield row.get_children()[0].get_children()

def main(sim = False):
	import numpy as np
	import matplotlib.pyplot as plt
	from subprocess import run
	import seaborn as sns
	import sys
	import scipy.interpolate as pol
	from scipy.optimize import curve_fit
	from colour import Color
	import sys
	import os.path
	import os



	L_ = np.array([40,60, 80,100])
	dT = 0.001
	Tstart = 2.25
	Tend = 2.3
	stable_logMCCs = 4
	stable_MCCs = int(10**stable_logMCCs)
	logMCCs = 4.5
	MCCs = int(10**logMCCs)
	TC = np.zeros(4)
	filenames = [f"paralell_L{L}_T0{Tstart}_T1{Tend}_dT{dT}_MCCs{stable_logMCCs}.dat" for L in L_]

	
	if sim:
		for L,filename in zip(L_,filenames):
			#run(f"rm ../data/{filename}".split())
			#run(f"../cpp/paralell.out {L} {MCCs} {stable_MCCs} {Tstart} {Tend+dT} {dT} {filename}".split())
			pass

	with sns.axes_style("darkgrid"):
		colors = ["red", "orange", "green", "blue"]
		fig, ax = plt.subplots(nrows=1, ncols=2, dpi=120)
		for i,(L,filename) in enumerate(zip(L_,filenames)):
			data = np.loadtxt(f"../data/{filename}")
			#E, M, E2, M2, Mabs, varE, varM, T, MCCs = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
			varE = data[:,5]
			T = data[:,7]
			sortIdx = np.argsort(T)
			
			T = np.array([T[idx] for idx in sortIdx])
			varE = np.array([varE[idx] for idx in sortIdx])/L**2
			CV = varE/T**2

			
			cs = pol.UnivariateSpline(T,CV,s=10000)
			Ts = np.linspace(T[0],T[-1],10000)

			print(Ts[np.argmax(cs(Ts))])
			TC[i] = Ts[np.argmax(cs(Ts))]
			ax[0].plot(Ts,cs(Ts) , color=Color(colors[i], luminance=0.4).get_rgb(), alpha=0.5,lw=2.2)
			ax[0].scatter(T, CV , color=colors[i], alpha=0.7,label=f"$L={L}$")
			ax[0].scatter(Ts[np.argmax(cs(Ts))], np.max(cs(Ts)), color=Color(colors[i], luminance=0.26).get_rgb(), marker="x", s=100, zorder=100)

		lines, labels = ax[0].get_legend_handles_labels()
		ax[0].legend(lines,labels,loc='upper left', ncol=2, fancybox=True, shadow=True,fontsize=10)
		fs = 14

		ax[0].set_ylabel("$C_V$",fontsize=fs)
		ax[0].set_xlabel("$T$ $[J/k]$",fontsize=fs)
		
		ax[0].set_xlim(2.25,2.3)
		ax[0].set_ylim(1.8,2.7)
		ax[0].set_title("$C_V(T)$ for $L = 40,60,80,100$.")



		L = np.array([40,60,80,100])

		linear_func = lambda x,a,b: a + b*x
		popt, pcov = curve_fit(linear_func, 1/L, TC)	
		m,b = np.polyfit(1/L, TC, 1)
		xf = np.linspace(1/L[-1], 1/L[0], 1000)
		print(b, pcov[1][1])




	

		colors = [Color(c, luminance=0.45).get_rgb() for c in ["red", "orange", "green", "blue"]]

		#fig, ax = plt.subplots(nrows=1, ncols=1, dpi=160)
		points = [0]*4
		for i,(L_,TC_,color) in enumerate(zip(L,TC,colors)):
			points[i] =ax[1].scatter(1/L_,TC_, marker="x", color=color, s=70, zorder=100)

		fit = ax[1].plot(xf, m*xf+b, c="k", lw=2)
		ax[1].set_xlim(0.0082,0.026)
		ax[1].set_ylim(2.276,2.29)
		ax[1].set_xlabel("$1/L$",fontsize=14)
		ax[1].set_ylabel("$T_C$ $[J/k]$",fontsize=14)
		ax[1].set_title("Linear fit of $T_C$ for $L = 40,60,80,100$.")
		ax[1].set_xticks([0.01, 0.015,0.02,0.025])
		#ax[1].set_yticks([2.280, 2.290])
		l = ax[1].legend([tuple(points),fit[0]], ["Datapoints, $T_C(L)$", "Fit, $T_C(L) = %.3f/L + %.3f$"%(m,b)],scatterpoints=4, loc="upper left", fancybox=True, shadow=True)


		handles_list = list(get_handle_lists(l))
		handles = handles_list[0] 
		for i in range(len(colors)):
			handles[i].set_edgecolors(colors[::-1])



		plt.show()
main(sim=False)