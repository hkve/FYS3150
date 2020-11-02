import numpy as np
import sys
import argparse
import subprocess

# A simulation can be run from terminal, see README for usage or run python3 master.py -h
parser = argparse.ArgumentParser(description="""
Solve a solar system simulation.
""")
parser.add_argument('dt', metavar='dt', type=float, help='Time step in log10(year) for simulation')
parser.add_argument('N', metavar='N', type=float, help='Number of integration points in log10')
parser.add_argument('-method', metavar='METHOD', type=str, help='Integration method to use. Must be "euler" or "verlet". Default: "verlet".', choices=["euler", "verlet"], default="verlet")
parser.add_argument('-beta', type=float, default=2.0,help='Change the inverse proportionality of gravity. Must be in range [2,3]. Default = 2')
parser.add_argument('-sys', metavar="file", default="sys.dat",help="Name of init file in initdData/ dir. Default: sys.dat")
parser.add_argument('-out', metavar="file", default="sys.out",help="Name of file in data/ dir where simulation results are stored. Default: sys.out")
parser.add_argument('-Nwrite', metavar='points', type=int, help='Numbers of data points to be stored/written. Defaults to N. NB, not logarithmic!', default=-1)
parser.add_argument('--GR', action='store_true', help='Do simulation with general relativity correction term.')
parser.add_argument('--fixSun', action='store_true', help='Do simulation with sun fixed at 0,0,0. (Sun must be first element in init file!)')
parser.add_argument('--compile', action='store_true', help='Compile main.cpp to "main.exe" before running')
parser.add_argument('--q', action='store_true', help='Run quietly')

def simulate(N, dt, method="verlet", beta=2, Nwrite=-1, GR=False, fixSun=False, compile=False, quiet=False, sys="sys.dat", out="sys.out"):
    """
    This function is called by all other python programs to run a simulation using main.exe
        Args:
            N: (int/float) log10 of the number of integration steps
            dt: (int/float) log10 of the step length [AU/yr]
            Nwrite: (int)  how many points to write
            GR: (bool) If GR term should be included in the calculations
            fixSun (bool) If the Sun should be fixed at the center
            complile: (bool) If main.cpp should be compiled
            quiet: (bool) If master.py should display output while preforming a simulation
            sys: (string) file containing initial conditions, located in initData dir
            out: (string) file to write data from simulation, located in data dir 
    """
    if not (2<= beta <= 3):
        if not quiet:
            print(f"ERROR; beta = {beta} is not in range [2,3]. Terminating.")
        sys.exit()
    if compile:
        if not quiet:
            print("Compiling...")
        subprocess.run(f"g++ -o main.exe main.cpp -O3".split())

    booldic = {True:1, False:0}
    q = booldic[quiet]

    if not quiet:
        if Nwrite  <= 0:
            NwriteStr = "N"
        else:
            NwriteStr = str(Nwrite)
        print(f"Solving for dt: {'%e' %10.0**dt}, N: {'%i' %10.0**N}, method: {method}, Nwrite: {NwriteStr}, beta: {beta}, GR:{GR}, hotel: Trivago, read: {sys}, write: {out}.")
    GR = booldic[GR]
    fixSun = booldic[fixSun]
    method = {"euler":0,"verlet":1}[method]

    if not (1< Nwrite <= 10**N):
        Nwrite = 10**N
  
    subprocess.run(f"./main.exe {sys} {out} {Nwrite} {'%e' %dt} {'%e' %N} {method} {beta} {GR} {fixSun} {q}".split())
    if not q:
        print("Done!")

if __name__ == "__main__":
    args = parser.parse_args()
    simulate(args.N, args.dt, args.method, args.beta, args.Nwrite, args.GR, args.fixSun, args.compile, args.q, args.sys, args.out)

