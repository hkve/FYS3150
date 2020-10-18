import numpy as np
import matplotlib.pyplot as plt 
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description="""
Solve a solar system simulation.
""")
parser.add_argument('dt', metavar='dt', type=float, help='Time step (in days) for simulation')
parser.add_argument('N', metavar='N', type=int, help='Number of integration points')
                
parser.add_argument('-method', metavar='METHOD', type=str, help='Integration method to use. Must be "euler" or "verlet". Default: "verlet".', choices=["euler", "verlet"], default="verlet")
parser.add_argument('-beta', type=float, default=2.0,help='Variable for 3e. Must be in range [2,3]. Default = 2')
parser.add_argument('-sys', metavar="file", default="sys.dat",help="Name of file where initial system is stored. Default: sys.dat")
parser.add_argument('-out', metavar="file", default="sys.out",help="Name of file where simulation results are stored. Default: sys.out")
parser.add_argument('-dpts', metavar='DPTS', type=int, help='Numbers of data points to be stored. Default: same as N', default=-1)
parser.add_argument('-time', metavar='TYPE', type=str, choices =["days", "years"] ,help='Decide if the time units of dt and init vel are in days or years. Default: days', default="days")
parser.add_argument('--GR', action='store_true', help='Do simulation with general relativity correction term.')
parser.add_argument('--compile', action='store_true', help='Compile main.cpp to "main.exe" before running')
parser.add_argument('--q', action='store_true', help='Run quiet')


args = parser.parse_args()





if __name__ == "__main__":
    
    if not (2<= args.beta <= 3):
        if not args.q:
            print(f"ERROR; beta = {args.beta} is not in range [2,3]. Terminating.")
        sys.exit()
    if args.compile:
        if not args.q:
            print("Compiling...")
        subprocess.run(f"g++ -o main.exe main.cpp -O3".split())
    dpts = args.dpts
    if not (0 < args.dpts <= args.N):
        dpts = args.N
    if not args.q:
        print(f"Solving for dt:{args.dt}, N:{args.N}, method:{args.method}, beta:{args.beta}, GR:{args.GR}")
        print(f"read: {args.sys}, write: {args.out}")
    GR = {True: 1, False: 0}[args.GR]
    q = {True: 1, False: 0}[args.q]
    method = {"euler":0,"verlet":1}[args.method]
    subprocess.run(f"./main.exe {args.sys} {args.out} {dpts} {args.dt} {args.N} {method} {args.beta} {GR} {args.time} {q}".split())
    if not args.q:
        print("Done!")
    