import numpy as np
import matplotlib.pyplot as plt
from subprocess import run
import time

from file_reader import read_data_file
from getInitialConditions import setInitialConditions
from master import simulate


def forward_simulation(dt, N, system_dict, GR=True):
    # forwards the simulation according to dt and N
    # only initial and final datapoints are stored
    setInitialConditions("precession_init.dat", system_dict)
    simulate(N, dt, Nwrite=2, GR=GR, sys="precession_init.dat", out="precession.dat", quiet=True)
    # run(f'python3.8 master.py {dt} {N} -sys initData/precession_init.dat -out precession.dat -Nwrite 2 -time years -method verlet --q{GR} --log'.split() )
    


def read_final_state(filename="precession.dat"):
    # reads the final datapoints of precession.dat and creates a system-dictionary of them
    system = read_data_file(filename)
    rM = system['Mercury'].r[:,-1]
    vM = system['Mercury'].v[:,-1]
    rS = system['Sun'].r[:,-1]
    vS = system['Sun'].v[:,-1]

    system_dict = {
        "Sun": [rS[0], rS[1], rS[2], vS[0], vS[1], vS[2]],
        "Mercury": [rM[0], rM[1], rM[2], vM[0], vM[1], vM[2]]
    }

    return system_dict


def getPerihelionAngle(dt, GR=True):
    # finds the perihelion angle of the last orbit of Mercury around the Sun that has been simulated
   
    GR_string = {True: "GR", False: "CLASSIC"}[GR]
    system = read_data_file(f"prec/final_{dt}_{GR_string}.dat")
    rM = system['Mercury'].r
    rS = system['Sun'].r
    rel_pos = rM-rS # relative position of Mercury and Sun
    rel_dist = np.linalg.norm(rel_pos, axis=0) # relative distance of Mercury and Sun

    peridx = rel_dist.argmin() # index value for the shortest distance between Mercury and Sun
    x_peri = rel_pos[0, peridx]
    y_peri = rel_pos[1, peridx]

    per_ang = np.rad2deg( np.arctan2(y_peri, x_peri) ) # perihelion angle of x-axis in degree

    return per_ang * 3600 # converts to arcseconds


def findPrecession(dt, GR=False, newsim=False):
    # Runs 100 one-year simulations adding up to 100 years forward, and finds the final perihelion angle.

    start = time.time() # timer

    # setting initial conditions for Mercury-Sun system with perihelion along x-axis
    system_dict = {"Sun": [0,0,0,0,0,0], "Mercury": [0.3075,0,0,0,12.44,0]}
    N =-dt # log10 of N, results in one-year simulation
    
    if newsim: # simulates from scratch
        for i in range(100):
            # simulates 1 years at a time
            forward_simulation(dt, N, system_dict, GR=GR)
            print(dt, i, GR) # just to keep track of simulations
            system_dict = read_final_state()

        # finally simulates just over one orbit, saving all the points to find the perihelion angle
        setInitialConditions("precession_orbit_init.dat", system_dict)
        N = np.log10(0.3) - dt # simulates 0.3 of a year, which is just over one orbit
        print("final")
        GR_string = {True:"GR", False:"CLASSIC"}[GR]
        simulate(N, dt, Nwrite=int(1e6), GR=GR, sys="precession_orbit_init.dat", out=f"prec/final_{dt}_{GR_string}.dat", quiet=True)

    # Gets the precession of the century simulated and prints
    precession_per_century = getPerihelionAngle(dt, GR=GR)

    print(f"The precession over one century (GR: {GR}) was {precession_per_century} arcseconds.")
    print(f"Took {time.time()-start:.1f} seconds to run")
    return precession_per_century



if __name__ == '__main__':
    dts = np.array([-8, -7, -6, -5,-4]) # array of stepsizes to check for

    classic = np.zeros(len(dts)) # empty array to save Newtonian precessions
    GR = np.zeros(len(dts)) # empty array to save GR precessions

    newsim = True # whether we are running new simulations or presenting already simulated data

    # plots the perihelion angles with and without the GR-correction term as a function of stepsize
    for i,dt in enumerate(dts):
        classic[i] = findPrecession(dt, GR=False, newsim=newsim)
        GR[i] = findPrecession(dt, GR=True, newsim=newsim)
    plt.plot(-dts, np.abs(GR), label="GR")
    plt.plot(-dts, np.abs(classic), label="classic")
    plt.grid()
    plt.legend()
    plt.show()
