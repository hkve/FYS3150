import numpy as np
import matplotlib.pyplot as plt
from subprocess import run

from colour import Color
import matplotlib as mpl
import seaborn as sns

from file_reader import read_data_file
from getInitialConditions import setInitialConditions
from master import simulate

def run_simulation(dt, T, v):
    # runs a single escape simulation and returns bool whether Earth escaped
    N = np.log10(T)-dt # log10 of N
    
    system_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,v,0]}

    setInitialConditions("escape_init.dat", system_dict)
    simulate(N=N, dt = dt, Nwrite=2, sys="escape_init.dat", out="escape.dat", fixSun=True, quiet=True)

    return check_escape()


def check_escape():
    # checks for a given simulation whether the Earth fulfills the escape condition
    system = read_data_file("escape.dat") # reads simulation data
    r = system["Earth"].r[:,-1] # grabs last position
    v = system["Earth"].v[:,-1] # grabs last veloicty

    v_rad = np.dot(r, v) / np.linalg.norm(r) # calculating the radial velocity

    escape = v_rad**2/np.dot(v, v) >= 0.99 # if the velocity is almost solely radial, Earth escaped
    return escape


def compute_escape_velocity(v0, dv, dt, T, max_iterations=1000):
    """
    For a given starting velocity and velocity increments and simulation parameters, this function
    simulates the Sun-Earth system and returns the escape velocity.
        Args:
            v0: (float) The first velocity to check for in AU/yr
            dv: (float) The incremental velocity increase after failed escape in AU/yr
            dt: (float) log10 of the time stepsize to run the simulation in yr.
            T: (float) The number of years to run each simulation for before checking the escape condition.
            max_iterations: (int) The number of attempts to escape before breaking the program and exiting.
    """
    # initializing values
    escape = False
    iterations = 0
    v = v0

    while not escape:
        iterations += 1
        if iterations > max_iterations:
            raise Exception(f"Earth was not able to escape within {max_iterations} iterations.")
        escape = run_simulation(dt, T, v)
        v += dv
    return v
    

def escapeVelocity(n=20, T = 100):
    """
    Code to plot different inital velocities in the v_y direction
        Args: 
            n: (int) Number of orbits/escapes
            T: (float/int) Total time to run over 
    """
    v = np.linspace(2*np.pi, 3*np.pi, n, endpoint=True)

    T = 100
    dt = -4    
    N = np.log10(T)-dt
    
    orbits = []
    for i in range(n):
        system_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,v[i],0]} # Store the different velocities
        setInitialConditions(f"escape_init_{i}.dat", system_dict)
        simulate(N=N, dt = dt, Nwrite=1000, sys=f"escape_init_{i}.dat", out=f"escape_{i}.dat", fixSun=True, quiet=True)
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
    v_escape = compute_escape_velocity(v0, dv, dt=-4, T=100)

    print(f"The escape velocity was ({v_escape:.2f} +- {dv:.2f}) AU/yr.")
    print(f"The theoretical escape velocity is {np.sqrt(8*np.pi**2):.2f} AU/yr")    
