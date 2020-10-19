import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions
import matplotlib.pyplot as plt
from subprocess import run


def forward_simulation(dt, N, system_dict, GR=" --GR"):
    # forwards the simulation according to dt and N
    # only initial and final datapoints are stored
    setInitialConditions("precession_init.dat", system_dict)
    run(f'python3.8 master.py {dt} {N} -sys initData/precession_init.dat -out precession.dat -Nwrite 2 -time years -method verlet --q{GR}'.split() )



def read_final_state():
    # reads the final datapoints of precession.dat and creates a system-dictionary of them
    system = read_data_file(f"precession.dat")
    rM = system['Mercury'].r[:,-1]
    vM = system['Mercury'].v[:,-1]
    rS = system['Sun'].r[:,-1]
    vS = system['Sun'].v[:,-1]

    system_dict = {"Sun": [rS[0], rS[1], rS[2], vS[0], vS[1], vS[2]], "Mercury": [rM[0], rM[1], rM[2], vM[0], vM[1], vM[2]]}

    return system_dict


def getPerihelionAngle(dt, system_dict, GR=" --GR"):
    # simulates approximately one orbit with dt, and finds the perihelion angle off the x-axis
    setInitialConditions("precession_orbit_init.dat", system_dict)

    N = int(0.3 / dt) # simulates 0.3 of a year, which is just over one orbit
    if N > 1e6:
        dpts = 1e6
    else:
        dpts = N
    run(f'python3.8 master.py {dt} {N} -sys initData/precession_orbit_init.dat -out precession_orbit.dat -Nwrite {dpts} -time years -method verlet --q{GR}'.split() )

    import time
    time.sleep(0.1)
    system = read_data_file("precession_orbit.dat")
    rM = system['Mercury'].r
    rS = system['Sun'].r
    rel_pos = rM-rS
    rel_dist = np.linalg.norm(rel_pos, axis=0)

    peridx = rel_dist.argmin()
    x_peri = rel_pos[0, peridx]
    y_peri = rel_pos[1, peridx]

    per_ang = np.rad2deg( np.arctan2(y_peri, x_peri) )

    return per_ang * 3600 # converts to arcseconds



if __name__ == '__main__':
    system_dict = {"Sun": [0,0,0,0,0,0], "Mercury": [0.3075,0,0,0,12.44,0]}
    dt = float(1e-6)
    N = int(1e7) # 1e7 is the largest array memory can handle

    for i in range(1):
        # simulates 10 years at a time
        forward_simulation(dt, N, system_dict, GR="")
        system_dict = read_final_state()

    precession_per_century = getPerihelionAngle(float(1e-6), system_dict, GR="")

    print(f"The precession over one century was {precession_per_century} arcseconds.")