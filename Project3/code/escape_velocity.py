import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions
import matplotlib.pyplot as plt
from subprocess import run

from colour import Color
import matplotlib as mpl
import seaborn as sns

def run_simulation(dt, T, v):
    N = int(T/dt + 1)

    system_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,v,0]}

    setInitialConditions("escape_init.dat", system_dict, fixedCoM = True)

    run(f'python master.py {dt} {N} -sys initData/escape_init.dat -out escape.dat -Nwrite 2 -time years -method verlet --q'.split() )

    return check_escape()


def check_escape():
    system = read_data_file("escape.dat")
    r = system["Earth"].r[:,-1]
    v = system["Earth"].v[:,-1]

    v_rad = np.dot(r, v) / np.linalg.norm(r)

    escape = v_rad**2/np.dot(v, v) >= 0.99 # if the velocity is almost solely radial, Earth escaped
    return escape


def compute_escape_velocity(v0, dv, dt, T, max_iterations=1000):
    escape = False
    iterations = 0
    v = v0
    while not escape:
        escape = run_simulation(dt, T, v)
        iterations += 1
        if iterations == max_iterations:
            raise Exception(f"Earth was not able to escape within {max_iterations} iterations.")
        v += dv
    return v
    

def plot_escape_velocity(n=20, T = 100):
    v = np.linspace(2*np.pi, 3*np.pi, n, endpoint=True)

    T = 100
    dt = 1e-4    
    N = int(T/dt + 1)

    orbits = []
    for i in range(n):
        """
        system_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,v[i],0]}
        setInitialConditions(f"escape_init_{i}.dat", system_dict, fixedCoM = False)
        run(f'python3 master.py {dt} {N} -sys initData/escape_init_{i}.dat -out escape_{i}.dat -Nwrite 1000 -time years -method verlet'.split())
        """
        system = read_data_file(f"escape_{i}.dat")
        orbits.append(system["Earth"].r)
        
    NoOfColors = n
    colors = list(Color("cyan").range_to(Color("orange"),NoOfColors)) 
    colors = [color.get_rgb() for color in colors]
    
    vticks = [round(v_i, 2) for v_i in v]
    
    fig, ax = plt.subplots(1,1)
    ax.set_facecolor('black')
    ax.set_xlabel("x [AU]", fontsize=15)
    ax.set_ylabel("y [AU]", fontsize=15)
    ax.set(xlim=(-50,2), ylim=(-26,26))
    for i, r in enumerate(orbits):
        ax.plot(r[0], r[1], color = colors[i], alpha=0.7)

    cmap = mpl.colors.ListedColormap(colors)
    norm = mpl.colors.Normalize(vmin = v[0],vmax = v[-1])
    cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical',ticks=vticks)
    cbar.ax.set_ylabel(r'   $v_y$', rotation=0, fontsize=17)
    cbar.ax.tick_params(labelsize=13)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.scatter(0,0, c="yellow", s=2)
    ax.text(0,15.5,'               $v_{esc}$', color="black", fontsize=17)
    ax.text(-49,6.5, "$v_{esc}=8.89$ [AU/yr]", color="white", fontsize=14)
    ax.arrow(-43,8.5,4,4,ec="white", fc="white", head_width=1)
    plt.show()

if __name__ == '__main__':
    v0 = 2*np.pi
    dv = 0.01
    v_escape = compute_escape_velocity(v0, dv, dt=1e-4, T=100)

    print(f"The escape velocity was ({v_escape:.2f} +- {dv:.2f}) AU/yr.")
    print(f"The theoretical escape velocity is {np.sqrt(8*np.pi**2):.2f} AU/yr")    
    