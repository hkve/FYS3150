import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import argparse

from stability import circularOrbit, ellipticalOrbits, benchmark, error
from jupiter_influence import sunEarthJupiter, radialDistance
from escape_velocity import escapeVelocity
from full_system import fullSystem
from modified_gravity import beta10yr, betaPosEnergy

mpl.use(plt.get_backend())
font = {'family' : 'normal',
        'size'   : 20}
mpl.rc('font', **font)


parser = argparse.ArgumentParser(description="""
Plot figures from the project
""")

# From stability
parser.add_argument('-circularOrbit', help='Plot stable circulat orbit of earth, with sun kept fixed',action='store_true',)
parser.add_argument('-ellipticalOrbits', help='Plot elliptical orbits of earth for different inital velocities, with sun kept fixed',action='store_true',)
parser.add_argument('-benchmark', help='Preforms a benchmarking of Euler and Velocity Verlet and plots the result',action='store_true',)
parser.add_argument('-error', help='Check energy conservation for Euler and Velocity Verlet and plots the result',action='store_true',)

# From jupiter_influence
parser.add_argument('-sunEarthJupiter', help='Plot Sun-Earth-Jupiter system with the mass of jupiter scaled by a facotr of 1000',action='store_true',)
parser.add_argument('-radialDistance', help='Plot the radial deviation from Earths orbit without Jupiter for different masses of jupiter',action='store_true',)

# From escape_velocity
parser.add_argument('-escapeVelocity', help='Plot a visual representation of escape velocity for different inital velocities',action='store_true',)

# From full_system
parser.add_argument('-fullSystem', help='Plot our own real system', action='store_true',)

# From modified_gravity
parser.add_argument('-beta10yr', help='Path of planet for beta=2.92 over 10 years',action='store_true',)
parser.add_argument('-betaPosEnergy', help='Plot path and energy of several beta-sytems', action='store_true',)


if __name__ == "__main__":
    args_ = parser.parse_args()
    args = []
    for arg in dir(args_)[::-1]:
        if arg == "_get_kwargs":
            break
        if eval(f"args_.{arg}"):
            args.append(arg)
    for arg in args:
        exec(f"{arg}()")