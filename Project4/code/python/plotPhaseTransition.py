def main(sim = False):
    import numpy as np
    import matplotlib.pyplot as plt
    from subprocess import run
    import seaborn as sns
    import sys
    L_ = np.array([40,60,80,100])
    #L_ = np.array([100])
    dT = 0.01
    Tstart = 2
    Tend = 2.3
    logMCCs = 5
    MCCs = int(10**logMCCs)

    filenames = [f"paralell_L{L}_T0{Tstart}_T1{Tend}_dT{dT}_MCCs{logMCCs}.dat" for L in L_]

    if sim:
        for L,filename in zip(L_,filenames):
            #run(f"rm ../data/{filename}".split())
            print(L)
            #run(f"../cpp/paralell.out {L} {MCCs} {Tstart} {Tend+dT} {dT} {filename}".split())


    with sns.axes_style("darkgrid"):
        colors = ["red", "orange", "green", "blue"]
        fig, axes = plt.subplots(nrows=2, ncols=2, dpi=100)
        for i,(L,filename) in enumerate(zip(L_,filenames)):
            data = np.loadtxt(f"../data/{filename}")
            E, M, E2, M2, Mabs, varE, varM, T, MCCs = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
            sortIdx = np.argsort(T)

            T = np.array([T[idx] for idx in sortIdx])
            E = np.array([E[idx] for idx in sortIdx])/L**2
            Mabs = np.array([Mabs[idx] for idx in sortIdx])/L**2
            varE = np.array([varE[idx] for idx in sortIdx])/L**2
            varM = np.array([varM[idx] for idx in sortIdx])/L**2
            CV = varE/T**2
            X = varM/T

            for ax, param in zip(np.array(axes).flatten(), [E,Mabs, CV, X]):
                ax.plot(T[:],param[:], "-o", color=colors[i], alpha=0.7,label=f"$L={L}$")
            print(colors[i], L)
            

           
        #axes[0,0].legend([f"$L={L}$" for L in L_])
        lines, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(lines,labels,loc='upper center', ncol=4, fancybox=True, shadow=True,fontsize=9)
        axes[0,0].set_title("$<E>$")
        axes[0,1].set_title("$< |M |>$")
        axes[1,0].set_title("$C_V$")
        axes[1,1].set_title("$\\chi$")
        plt.show()
main(sim=False)