import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions
import matplotlib.pyplot as plt
from subprocess import run
import time


def forward_simulation(dt, N, system_dict, GR=" --GR"):
    # forwards the simulation according to dt and N
    # only initial and final datapoints are stored
    setInitialConditions("precession_init.dat", system_dict)
    run(f'python3.8 master.py {dt} {N} -sys initData/precession_init.dat -out precession.dat -Nwrite 2 -time years -method verlet --q{GR} --log'.split() )
    


def read_final_state(filename="precession.dat"):
    # reads the final datapoints of precession.dat and creates a system-dictionary of them
    system = read_data_file(filename)
    rM = system['Mercury'].r[:,-1]
    vM = system['Mercury'].v[:,-1]
    rS = system['Sun'].r[:,-1]
    vS = system['Sun'].v[:,-1]

    system_dict = {"Sun": [rS[0], rS[1], rS[2], vS[0], vS[1], vS[2]], "Mercury": [rM[0], rM[1], rM[2], vM[0], vM[1], vM[2]]}

    return system_dict


def getPerihelionAngle(dt, GR=" --GR"):
    # simulates approximately one orbit with dt, and finds the perihelion angle off the x-axis
   
   
    system = read_data_file(f"prec/final_{dt}_{GR}.dat")
    #run(f'cp data/precession_orbit.dat data/prec/{GR}.{dt}.dat'.split())
    rM = system['Mercury'].r
    rS = system['Sun'].r
    rel_pos = rM-rS
    rel_dist = np.linalg.norm(rel_pos, axis=0)

    peridx = rel_dist.argmin()
    x_peri = rel_pos[0, peridx]
    y_peri = rel_pos[1, peridx]

    per_ang = np.rad2deg( np.arctan2(y_peri, x_peri) )

    return per_ang * 3600 # converts to arcseconds


def findPrecession(dt, GR=False, newsim=False):
    start = time.time()
    system_dict = {"Sun": [0,0,0,0,0,0], "Mercury": [0.3075,0,0,0,12.44,0]}
    N =-dt # 1e7 is the largest array memory can handle
    if GR:
        filename = f'GR.{dt}.dat'
    else:
        filename = f'CLASSIC.{dt}.dat'
    
    GR_s = {True:" --GR", False:""}[GR]
    if newsim:
        for i in range(100):
            # simulates 1 years at a time
            forward_simulation(dt, N, system_dict, GR=GR_s)
            print(dt, i, GR)
            system_dict = read_final_state()
        GR_f = {True:"GR", False:"CLASSIC"}[GR]
        run(f'cp data/precession.dat data/prec/99yr_{dt}_{GR_f}.dat'.split())
        system_dict = read_final_state(f"prec/99yr_{dt}_{GR_f}.dat")

        setInitialConditions("precession_orbit_init.dat", system_dict)
        N = np.log10(0.3 / 10**float(dt)) # simulates 0.3 of a year, which is just over one orbit
        print("final")
        run(f'python3.8 master.py {dt} {N} -sys initData/precession_orbit_init.dat -out prec/final_{dt}_{GR_f}.dat -time years -Nwrite {int(round(10**float(N)))} -method verlet --q{GR_s} --log'.split() )

    #else:
    #    system_dict = read_final_state(f"prec/final_{dt}_{GR_}.dat")
    GR_ = {True:"GR", False:"CLASSIC"}[GR]

    precession_per_century = getPerihelionAngle(dt, GR=GR_)

    print(f"The precession over one century (GR: {GR}) was {precession_per_century} arcseconds.")
    print(f"Took {time.time()-start:.1f} seconds to run")
    return precession_per_century
if __name__ == '__main__':
    dts = np.array([-8,-7,-6,-5,-4])#np.linspace(-5,-5,1, dtype=np.int16)[::-1]
    print(dts)
    classic = np.zeros(len(dts))
    GR = np.zeros(len(dts))

    newsim = False

    for i,dt in enumerate(dts):
        classic[i] = findPrecession(dt, GR=False, newsim=newsim)
        GR[i] = findPrecession(dt, GR=True, newsim=newsim)
    plt.plot(-dts, np.abs(GR), label="GR")
    plt.plot(-dts, np.abs(classic), label="classic")
    plt.grid()
    plt.legend()
    plt.show()
