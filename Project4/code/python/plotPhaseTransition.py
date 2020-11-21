def getkwarg(key, kwargs):
    try:
        if key in kwargs.keys():
            return kwargs[key]
    except:
        return None

def main(sim = False, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from subprocess import run
    import seaborn as sns
    import sys
    import scipy.interpolate as pol

    L_ = np.array([2,40])#40,60])#,80,100])
    L_ = [100]
    #L_ = np.array([100])
    dT = 0.05
    Tstart = 1.5
    Tend = 2.3
    stable_logMCCs = 4.3
    stable_MCCs = int(10**stable_logMCCs)
    logMCCs = 5
    MCCs = int(10**logMCCs)

    filenames = [f"paralell_L{L}_T0{Tstart}_T1{Tend}_dT{dT}_MCCs{stable_logMCCs}.dat" for L in L_]
    
    if sim:
        for L,filename in zip(L_,filenames):
            run(f"rm ../data/{filename}".split())
            run(f"../cpp/paralell.out {L} {MCCs} {stable_MCCs} {Tstart} {Tend+dT} {dT} {filename}".split())


    with sns.axes_style("darkgrid"):
        colors = ["red", "orange", "green", "blue"]
        fig, axes = plt.subplots(nrows=2, ncols=2, dpi=100)
        for i,(L,filename) in enumerate(zip(L_,filenames)):
            data = np.loadtxt(f"../data/{filename}")
            E, M, E2, M2, Mabs, varE, varM, T, MCCs = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
            sortIdx = np.argsort(T)

            T = np.array([T[idx] for idx in sortIdx])-dT
            E = np.array([E[idx] for idx in sortIdx])/L**2
            M = np.array([M[idx] for idx in sortIdx])/L**2
            Mabs = np.array([Mabs[idx] for idx in sortIdx])/L**2
            varE = np.array([varE[idx] for idx in sortIdx])/L**2
            varM = np.array([varM[idx] for idx in sortIdx])/L**2
            CV = varE/T**2
            X = varM/T


            for ax, param in zip(np.array(axes).flatten(), [E,Mabs, CV, X]):
                pops = 28
                pops = [1.5,1.540,1.620,1.63,1.64,1.66,1.68,1.74,1.76,1.83,1.94]
                T_ = list(T)[:]
                param_ = list(param)[:]
                for val in pops:
                    midx = T_.index(val-dT)
                    T_.pop(midx)
                    param_.pop(midx)

                cs = pol.CubicSpline(T_,param_)
                Ts = np.linspace(T_[0],T_[-1],1000)
                ax.plot(Ts,cs(Ts) , color=colors[i], alpha=0.7,label=f"$L={L}$")
                ax.plot(T_, param_ , "-o", color=colors[i], alpha=0.7,label=f"$L={L}$")

           
        #axes[0,0].legend([f"$L={L}$" for L in L_])
        lines, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(lines,labels,loc='upper center', ncol=4, fancybox=True, shadow=True,fontsize=9)
        axes[0,0].set_title("$<E>$")
        axes[0,1].set_title("$< |M |>$")
        axes[1,0].set_title("$C_V$")
        axes[1,1].set_title("$\\chi$")
        plt.show()
main(sim=False)