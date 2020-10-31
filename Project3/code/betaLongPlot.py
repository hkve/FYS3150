import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions as sic
import matplotlib.pyplot as plt
import matplotlib as mpl
from subprocess import run
mpl.use(plt.get_backend())
font = {'family' : 'normal',
        'size'   : 20}
mpl.rc('font', **font)


simYears = 10 # how many years to simulate

dt = -7
N = np.log10(simYears/10**dt)
Nwrite = int(10**N/1000)
beta = 2.92


body_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,5,0]}
#sic("varyingBeta.dat", body_dict, fixedCoM=True)
 
#run(f'python3.8 master.py {_dt} {_N} -sys initData/varyingBeta.dat -out beta_{beta}.out -time years -Nwrite {Nwrite} -beta {beta} --log'.split() )


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
