import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions as sic
import matplotlib.pyplot as plt
from subprocess import run
import matplotlib as mpl
from colour import Color
from subprocess import run
import pandas as pd
mpl.use(plt.get_backend())
font = {'family' : 'normal',
        'size'   : 20}

mpl.rc('font', **font)

G_ = 6.67408e-11
AU = 1.496e+11
ME = 5.972e+24
G = G_/AU**3*(3600*24*365)**2*ME
c = 299792458/AU*3600*24*365

body_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,2*np.pi,0]}
dt = np.arange(-7,-1)[::-1]
print(dt)
N = -dt+1
print(N)

sic("stability.dat", body_dict)
for dt_, N_ in zip(dt,N):
    pass
    run(f'python3.8 master.py {dt_} {N_} -sys initData/stability.dat -out {-dt_}.{N_}.out -Nwrite {int(1e5)} -time years --log'.split() )

colors = list(Color("cyan").range_to(Color("orange"),len(dt))) 
colors = [str(color) for color in colors]
with plt.style.context("seaborn-darkgrid"):
    f, ax = plt.subplots(1,1, dpi=100,frameon=True)

   
    dic = {r"$\beta$":[r"$\log_{10} \sigma_{\mathcal{E}}$",r"$\log_{10}\sigma_{\ell}$"]}
    for idx,(dt_,N_) in enumerate(zip(dt,N)):
        if idx == 0:
            continue
        system = read_data_file(f"{-dt_}.{N_}.out")
        r1 = system["Earth"].r
        r0 = system["Sun"].r
        v1 = system["Earth"].v
        v0 = system["Sun"].v
        dist = np.linalg.norm(r1-r0, axis=-2)
        
        
        ax.plot(np.linspace(dt[0]*N[0], dt[-1]*N[-1], len(dist)), dist,color=colors[idx], lw=2, alpha=0.4)
  
    ax.set_facecolor((0.15, 0.15, 0.15))

    ax.set_xlabel("$x$ [AU]")
    ax.set_ylabel("$y$ [AU]")
    #ax.set_ylim(top=1.1, bottom=-1.1)
    #ax.set_xlim(-1.1, 1.1)
    cmap = mpl.colors.ListedColormap(colors)
    norm = mpl.colors.Normalize(vmin=dt[0],vmax=dt[-1])
    cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical')# ticks=[int(a) for a in np.linspace(betas[0], betas[-1],100)] )
    cbar.ax.set_ylabel(r'   $\beta$', rotation=0, fontsize=20)
    legend = ax.legend(fancybox=True, framealpha=1, shadow=True, borderpad=0.4, frameon = True,loc='lower right')
    legend.get_frame().set_facecolor('grey')
    ax.set_title(r"Earth-Sun system for $\beta \in [2,3]$")
plt.grid()
#plt.xticks([-1,0,1])
#plt.yticks([-1,0,1])
plt.legend()
plt.show()
