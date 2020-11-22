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
    threads = 4
    filenames = [f"paralell_L{L}_T0{Tstart}_T1{Tend}_dT{dT}_MCCs{stable_logMCCs}.dat" for L in L_]

    
    if sim:
        for L,filename in zip(L_,filenames):
            #run(f"rm ../data/{filename}".split())
            #run(f"../cpp/paralell.out {L} {MCCs} {stable_MCCs} {Tstart} {Tend} {dT} {threads} {filename}".split())
            pass

    with sns.axes_style("darkgrid"):
        colors = ["red", "orange", "green", "blue"]
        fig, axes = plt.subplots(nrows=2, ncols=2, dpi=120)
        for i,(L,filename) in enumerate(zip(L_,filenames)):
            data = np.loadtxt(f"../data/{filename}")
            E, M, E2, M2, Mabs, varE, varM, T, MCCs = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
            sortIdx = np.argsort(T)

            T = np.array([T[idx] for idx in sortIdx])
            E = np.array([E[idx] for idx in sortIdx])/L**2
            M = np.array([M[idx] for idx in sortIdx])/L**2
            Mabs = np.array([Mabs[idx] for idx in sortIdx])/L**2
            varE = np.array([varE[idx] for idx in sortIdx])/L**2
            varM = np.array([varM[idx] for idx in sortIdx])/L**2
            CV = varE/T**2
            X = varM/T


            for j,(ax, param) in enumerate(zip(np.array(axes).flatten(), [E,Mabs, CV,X])):
                
                cs = pol.UnivariateSpline(T_,param_,s=10000)
                Ts = np.linspace(T_[0],T_[-1],10000)
               
                ax.plot(Ts,cs(Ts) , color=Color(colors[i], luminance=0.4).get_rgb(), alpha=0.7,lw=2.2)
                ax.scatter(T_, param_ , color=colors[i], alpha=0.7,label=f"$L={L}$",)

        #axes[0,0].legend([f"$L={L}$" for L in L_])
        lines, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(lines,labels,loc='upper center', ncol=4, fancybox=True, shadow=True,fontsize=12)
        fs = 14
        axes[0,0].set_ylabel("$\\langle E\\rangle/L^2$",fontsize=fs)
        axes[0,1].set_ylabel("$\\langle |M |\\rangle/L^2$",fontsize=fs)
        axes[1,0].set_ylabel("$C_V/L^2$",fontsize=fs)
        axes[1,1].set_ylabel("$\\chi/L^2$",fontsize=fs)
        axes[1,1].set_xlabel("$T$ $[J/k]$",fontsize=fs)
        axes[1,0].set_xlabel("$T$ $[J/k]$",fontsize=fs)
        for i in range(2):
            for j in range(2):
                axes[i,j].set_xlim(2.25,2.3)
        axes[0,0].set_ylim(-1.47,-1.34)
        axes[1,0].set_ylim(1.8,2.6)
        axes[1,1].set_ylim(-5, 173)
        plt.show()
main(sim=False)