import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions, getInitialCondition
import matplotlib.pyplot as plt
from subprocess import run
from colour import Color
import matplotlib as mpl
mpl.use(plt.get_backend())
font = {'family' : 'normal',
        'size'   : 18}
mpl.rc('font', **font)



N = 7
dt = np.log10(10**(-7)*25)
Nwrite = int(1e4)

bodynames = ["Sun", "Mercury", "Venus","Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"]
bodycolors = ["yellow", "#86511D", "#F7C3F4", "#0EEB58", "red", "orange", "#FCDB0A", "aqua", "blue", "grey"]
#getInitialCondition("fullsystem.dat", bodies=bodynames, fixedCoM=False, timeFormat = "years", date="2018-03-14")
#run(f"python3.8 master.py {dt} {N} --log -time years -Nwrite {Nwrite} -sys initData/fullsystem.dat -out fullsystem.out --GR".split())

system = read_data_file("fullsystem.out")

with plt.style.context("seaborn-darkgrid"):
    fig = plt.figure(figsize=(13,8), dpi=130)
    ax = fig.gca(projection='3d')

    for i,bodyname in enumerate(bodynames):
        body = system[bodyname]
        r = body.r
        v = body.v
        
        h = np.cross(r, v,axis=-2)
        inc = np.mean(np.rad2deg(np.arccos(h[-1,:]/np.linalg.norm(h, axis=-2)))[-1]) 
        print(f"Inclination {bodyname}: {round(inc,3)} degr")
        x,y,z = body.r[0,:], body.r[1,:], body.r[2,:]
        
        ax.plot(x,y,z, lw=1, color = bodycolors[i], alpha=0.65,  zorder=-1)
        ax.scatter(x[-1], y[-1], z[-1], color = bodycolors[i], zorder=10,s=18)
        ax.text(x[-1], y[-1], z[-1], bodyname, color=Color(bodycolors[i], luminance=0.95).rgb , fontsize=12, zorder=1000)
    
 
    ax.set_axis_off()
    ax.set_facecolor("black")

    ax.set_title("All bodies over 25 years", color="white")
    lim = 20 # plot dimensions (AU)
    zscale = 20 # scaling of z-axis to show inclination
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_zlim(-lim/zscale,lim/zscale)
    plt.show()


