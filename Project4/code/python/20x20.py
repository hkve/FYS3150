import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from subprocess import run
import sys
L = 20
MCCs_start = 10
MCCs_end = 10000
dMCCs = 100
T_ = [1,2.4]
initSpin_ = [0,1]
sim = not True
MCCs = np.array([1,5,10,50,100,200,300,400,500,600,700,800,900,1000,1200,1500,1800,2200,2600,3000,3400,5000,6000,7000,8000,9000,10000,20000,30000,40000,80000,100000])
E_ = []
Mabs_ = []

for T in T_:
    for initSpin in initSpin_:
        filename = f"../data/20x20_{T}_{initSpin}.out"
        if sim:
            try:
                run(f"rm {filename}".split())
            except:
                print(f"failed to remove {filename}")
            for MCC in MCCs:
                print(T, initSpin, MCC)
                run(f"../cpp/20x20_lattice.out {T} {MCC} {initSpin} {filename}".split())
        data = np.loadtxt(filename)
        E, M, E2, M2, Mabs, varE, varM, T__ = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7]
        E_.append(E.copy())
        Mabs_.append(Mabs.copy())

spinLabel = ["Random", "Ordered"]
Tcolors = ["blue", "red"]
spinMarker = ["solid", (0,(5,1))]

with sns.axes_style("darkgrid"):
    fig = plt.figure(dpi=200)
    axes = fig.subplots(nrows=1,ncols=2, sharex=True)
    foo = 0
    for i,T in enumerate(T_):
        for j,initSpin in enumerate(initSpin_):
            axes[0].plot(MCCs, E_[foo]/L**2, linestyle=spinMarker[initSpin] ,lw=3, color = Tcolors[i], markersize=4, alpha=0.5)#,label=f"$T = {T}$, {spinLabel[initSpin]} spins")
            axes[0].plot(0, linestyle=spinMarker[initSpin] ,lw=1.7, color = Tcolors[i],label=f"$T = {T}$, {spinLabel[initSpin]} spins")

            axes[1].plot(MCCs, Mabs_[foo]/L**2, linestyle=spinMarker[initSpin], lw=3, color=Tcolors[i], alpha=0.5)#, label=f"$T = {T}$, {spinLabel[initSpin]} spins")
            foo += 1
    axes[0].set_ylim(-2.1, -0.5)
    axes[1].set_ylim(-0.1,1.1)
    lines, labels = fig.axes[0].get_legend_handles_labels()
    fig.legend(lines,labels,loc='upper center', ncol=2, fancybox=True, shadow=True,fontsize=9)#, bbox_to_anchor=(0.5, 1.16),ncol=2, fancybox=True, shadow=True, fontsize=9)
    axes[1].set_xlabel("MCCs", fontsize=12)
    axes[0].set_xlabel("MCCs", fontsize=12)

    axes[0].set_ylabel("$\\langle E \\rangle / L^2$", fontsize=13)
    axes[1].set_ylabel("$ \\left <|M|\\right > / L^2 $", fontsize=13)
plt.show()