import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions as sic
from getInitialConditions import getInitialCondition as gic
import matplotlib.pyplot as plt
import matplotlib as mpl
from master import simulate
from colour import Color
import pandas as pd

def beta10yr():
    simYears = 10 # how many years to simulate

    dt = -7
    N = np.log10(simYears/10**dt)
    Nwrite = 10000
    beta = 2.92


    body_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,5,0]}
    sic("varyingBeta.dat", body_dict, fixedCoM=True)

    simulate(N=N, dt = dt, beta=beta, Nwrite=Nwrite, sys="varyingBeta.dat", out=f"beta_{beta}.out")

    with plt.style.context("seaborn-darkgrid"):
        f, ax = plt.subplots(1,1, dpi=100,frameon=True)
        
    
        system = read_data_file(f"beta_{beta}.out")
        x = system["Earth"].r[0,:]
        y = system["Earth"].r[1,:]
        cmap=plt.get_cmap('jet')
        points = np.array([x,y]).T.reshape(-1, 1, 2)

        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = mpl.collections.LineCollection(segments, array=np.linspace(0.0, 1.0, len(x)),  cmap=cmap, norm=plt.Normalize(0.0, 1.0), linewidth=1)

    
        ax.add_collection(lc)
        cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=plt.Normalize(1.0, 10.0), cmap=cmap), ax=ax, orientation='vertical', ticks=list(range(1,11)))
        cbar.ax.set_ylabel(r'       year', rotation=-36, fontsize=20)

        ax.set_facecolor((0.15, 0.15, 0.15))

        ax.set_xlabel("$x$ [AU]")
        ax.set_ylabel("$y$ [AU]")
        ax.set_ylim(top=1.1, bottom=-1.1)
        ax.set_xlim(-1.1, 1.1)
        
        ax.set_title(r"Earth-Sun system for $\beta = $" + f"${beta}$")
        plt.grid()
        plt.xticks([-1,0,1])
        plt.yticks([-1,0,1])
        plt.show()

def betaPosEnergy():
    GM0 = 4*np.pi**2
    N = 7
    dt = np.log10(10**(-7)*0.5)
    bnums = 10
    betas = np.linspace(2,3.0, bnums, endpoint=False)

    body_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,5,0]}

    sic("varyingBeta.dat", body_dict, fixedCoM=True)
    for beta in betas:
        simulate(N=N, dt=dt, sys="varyingBeta.dat", out=f"{round(beta,3)}.out", Nwrite=int(1e4), beta=beta)
        pass
    colors = list(Color("cyan", luminance=0.5).range_to(Color("orange", luminance=0.5),bnums)) 
    colors = [color.get_rgb() for color in colors]
    colors_dark = list(Color("cyan", luminance=0.45).range_to(Color("orange", luminance=0.45),bnums)) 
    colors_dark = [color.get_rgb() for color in colors_dark]
    with plt.style.context("seaborn-darkgrid"):
        fr, axr = plt.subplots(1,1, dpi=100,frameon=True)
        fE, axE = plt.subplots(1,1, dpi=100,frameon=True)
        ls, Es, Eadjs = [], [], []
        stdl = np.zeros(len(betas))
        stdE = np.zeros(len(betas))
        stdEadj = np.zeros(len(betas))


        dic = {r"$\beta$":[], r"$\sigma_{\ell}$":[],  r"$\bar{\ell}$":[], r"$\sigma_{\mathcal{E}}$":[], r"$\bar{\mathcal{E}}$":[]}

        for idx,beta in enumerate(betas):
            beta_ = round(beta,3)
            system = read_data_file(f"{beta_}.out")
            r1 = system["Earth"].r
            r0 = system["Sun"].r
            v1 = system["Earth"].v
            v0 = system["Sun"].v

            l = np.linalg.norm(np.cross(r1-r0, v1-v0,axis=-2), axis=-2)
            ls.append(l)
            stdl[idx] = np.std(l)
            
            E = 1/2*np.linalg.norm(v1-v0,axis=-2)**2 -1/(beta - 1)*GM0*1/((np.linalg.norm(r1-r0, axis=-2))**(beta-1))
            Es.append(E)
            stdE[idx] = np.std(E)
        
            maxidx = np.argmax(E)
            Eadj = np.append(E[:maxidx-1000], E[maxidx+1000:])
            Eadjs.append(Eadj)
            stdEadj = np.std(Eadj)
            
            dic[r"$\beta$"].append(beta_)
            dic[r"$\sigma_{\ell}$"].append("%.2e" %np.std(l))
            dic[r"$\bar{\ell}$"].append("%.2f" %np.mean(l))
            dic[r"$\sigma_{\mathcal{E}}$"].append("%.2e" %np.std(E))
            dic[r"$\bar{\mathcal{E}}$"].append("%.2f" %np.mean(E))
            
            
            axE.plot(np.linspace(0,0.5, len(E)),E, color=colors_dark[idx], lw=3, zorder=1000 -idx)#, label=r"$\beta = $"+f"${beta_}$" )
            
            axr.plot(np.array(system["Earth"].r[0,:]), np.array(system["Earth"].r[1,:]), color=colors[bnums-idx-1], lw=2, alpha=0.7)

        
        df = pd.DataFrame(dic)
        print(df.to_latex(index=False, escape=False, label="std"))

        cmap = mpl.colors.ListedColormap(colors)
        norm = mpl.colors.Normalize(vmin=betas[0],vmax=betas[-1])
        cbar = axr.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axr, orientation='vertical',ticks=betas )
        cbar.ax.set_ylabel(r'   $\beta$', rotation=0, fontsize=20)
        cbar = axE.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axE, orientation='vertical',ticks=betas )
        cbar.ax.set_ylabel(r'   $\beta$', rotation=0, fontsize=20)

        axr.set_xlabel("$x$ [AU]")
        axr.set_ylabel("$y$ [AU]")
        axr.set_ylim(top=1.1, bottom=-1.1)
        axr.set_xlim(left=-1.1, right=1.1)
        axr.grid(True)
        axr.set_xticks([-1,0,1])
        axr.set_yticks([-1,0,1])
        axr.set_title(r"Earth-Sun system for $\beta \in [2,3]$")
        axr.set_facecolor((0.15, 0.15, 0.15))

    
        axE.set_yticks(ticks=[-30, -20, -10, 0 ])
        axE.set_xticks(list(np.linspace(0,0.5, 3, endpoint=True)))
        axE.set_title("Energy in Earth-Sun system over 1/2 year")
        axE.set_xlabel("time, years")
        axE.set_ylabel(r"$\mathcal{E}$", rotation=0, fontsize=20)
        axE.set_ylim(-31, 1)
        plt.show()