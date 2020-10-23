import numpy as np
import matplotlib.pyplot as plt
from file_reader import read_data_file
from getInitialConditions import getInitialCondition
from subprocess import run


# solar_system = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
solar_system = ["Sun", "Mercury", "Earth"]
getInitialCondition("FSS_init.dat", bodies=solar_system, fixedCoM=True)


T = 365.25
dt = float(1e-4)
N = int(T/dt)
run(f'python master.py {dt} {N} -sys initData/FSS_init.dat -out FSS.dat -Nwrite 10000 -time days -method verlet'.split() )

system = read_data_file("FSS.dat")


fig, ax = plt.subplots()

NonBodies = ["dt", "N", "N_write", "method", "time"]
for body in system:
    if body not in NonBodies:
        x = system[body].r[0]
        y = system[body].r[1]
        ax.plot(x, y, label=body)

ax.axis('equal')
ax.legend()
plt.show()