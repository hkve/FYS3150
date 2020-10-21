import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions
import matplotlib.pyplot as plt
from subprocess import run



def run_simulation(dt, T, v):
    N = int(T/dt + 1)

    system_dict = {"Sun": [0,0,0,0,0,0], "Earth": [1,0,0,0,v,0]}

    setInitialConditions("escape_init.dat", system_dict)

    run(f'python master.py {dt} {N} -sys initData/escape_init.dat -out escape.dat -Nwrite 2 -time years -method verlet --q'.split() )

    return check_escape()


def check_escape():
    system = read_data_file("escape.dat")
    r = system["Earth"].r[:,-1]
    v = system["Earth"].v[:,-1]

    v_rad = np.dot(r, v) / np.linalg.norm(r)

    escape = v_rad**2/np.dot(v, v) >= 0.99 # if the radial velocity is greater than the angular, Earth escaped
    return escape


def compute_escape_velocity(v0, dv, dt, T, max_iterations=1000):
    escape = False
    iterations = 0
    v = v0
    while not escape:
        escape = run_simulation(dt, T, v)
        iterations += 1
        if iterations == max_iterations:
            raise Exception("Earth was not able to escape for any of the teste velocities.")
        v += dv
    return v
    



if __name__ == '__main__':
    v0 = 2*np.pi
    dv = 0.01
    v_escape = compute_escape_velocity(v0, dv, dt=1e-4, T=100)

    print(f"The escape velocity was ({v_escape:.2f} +- {dv:.2f}) AU/yr.")
    print(f"The theoretical escape velocity is {np.sqrt(8*np.pi**2):.2f} AU/yr")