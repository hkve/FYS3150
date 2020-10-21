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


def compute_escape_velocity(v_array, dt, T):
    escape = False
    for v in v_array:
        escape = run_simulation(dt, T, v)
        if escape:
            return v
    raise Exception("Earth was not able to escape for any of the teste velocities.")



if __name__ == '__main__':
    v_array = np.linspace(5, 15, 250)
    dv = v_array[1] - v_array[0]
    v_escape = compute_escape_velocity(v_array, dt=1e-4, T=100)

    print(f"The escape velocity was ({v_escape:.2f} +- {dv:.2f}) AU/yr.")
    print(f"The theoretical escape velocity is {np.sqrt(8*np.pi**2):.2f} AU/yr")