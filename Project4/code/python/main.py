import numpy as np
import sys, argparse, subprocess, argcomplete, os


import plot20x20, plotExpValues2x2, plotMCS2x2, plotAcceptedFlips, plotEnergyHistogram


class fakeStr(str):
  def __init__(self, val):
    self.val = val
  def __str__(self):
    return ""

class CapitalisedHelpFormatter(argparse.HelpFormatter):
  def add_usage(self, usage, actions, groups, prefix=None):
    return None
    
class FakeArgumentParser(argparse.ArgumentParser):
  def error(self, message):
    self.print_help(sys.stderr)

parser = FakeArgumentParser(description="Usage: main.py {prog} [args]", add_help=False, formatter_class=CapitalisedHelpFormatter)
parser._positionals.title = 'prog'
parser._optionals.title = 'args'

parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
argparse._HelpAction(option_strings=['-h', '--help'], dest='help', default='==SUPPRESS==', help='Show this help message and exit.')

parser.add_argument('type',type=str,  choices = ["2x2", "20x20", "parallel", "plot-exp-2x2", "plot-EM-20x20", "plot-MCCs-2x2", "plot-accepted-flips", "plot-E-hist"],help="Follow program name by optional program specific parameters listed below")
parser.add_argument("-MCCs", metavar="",type=int, help="int. Number of MCCs to perform")
parser.add_argument("-T", metavar="",type=float, help="float. System temperature")
parser.add_argument("-T0", metavar = "", type=float, help="float. Start temperature")
parser.add_argument("-T1", metavar = "", type=float, help="float. End temperature")
parser.add_argument("-dT", metavar = "", type=float, help="float. Temperature step")
parser.add_argument("-initSpin", metavar="", type=int, help="int. Init spin configuration. 1: all up, -1: all down, 0: random")
parser.add_argument("-filename", metavar="", type=str, help="str. Output filename")
parser.add_argument("-L", metavar="", type=str, help="Grid size")
parser.add_argument("-saveAfter", metavar="", type=str, help="When the grid energy becomes stable")
parser.add_argument("-threads", metavar="", type=str, help="Number of threads to use if parallel")
parser.add_argument("--sim", action="store_true", help="Re-run simulations in plotting code with default parameters. May take time")
parser.add_argument("--compile", action="store_true", help="Compile corresponding c++ program")

comp = parser.add_argument_group("If you simply want to compile all programs, write")
comp.add_argument(fakeStr("-fakecompile"), metavar = "compile", required=False, action="append", help="Compiles all C++ programs")

progparams = parser.add_argument_group("Program parameters (specifics) \nFrom here you can run the C++ programs")
progparams.add_argument(fakeStr("-fake2x2"), metavar = "2x2 [-T0] [-T1] [-dT] [-MCCs] [-filename] [--compile]", required=False, action="append", help="2x2 lattice from temperature T0 to T1 with step dT. Defaults: T0=1, T1=4, dT=0.1, MCCs=100000, filename=2x2EXP.out")
progparams.add_argument(fakeStr("-fake20x20"), metavar = "20x20 [-T] [-MCCs] [-initSpin] [-filename] [--compile]", required=False, action="append", help="20x20 lattice at temp. T. Defaults: T=1, MCCs=100000, initSpin = 0, filename=20x20.out")
progparams.add_argument(fakeStr("-parallel"), metavar="parallel [-L] [-MCCs] [-saveAfter] [-T0] [-T1] [-dT] [threads] [-filename] [--compile]", required=False, action="append", help="Runs different grid sizes in temperature range [T0, T1] with dT spacing. Defaults L=20 MCCs=1000000 saveAfter=10**4.2, T0=2, T1=2.3, dT=0.1, threads=2, filename=parallell_run.dat")

plots = parser.add_argument_group("Program plots \nFrom here you can recreate the plots as seen in the report")
plots.add_argument(fakeStr("-fakeplot2x2"), metavar = "plot-exp-2x2 [--sim] [--compile]", required=False, action="append", help="Plots exp. values for a 2x2 lattice")
plots.add_argument(fakeStr("-fakeplot2x2_MCCs"), metavar = "plot-MCCs-2x2 [--sim] [--compile]", required=False, action="append", help="Plots exp. values vs MCCs for a 2x2 lattice")
plots.add_argument(fakeStr("-fakeplot20x20"), metavar = "plot-EM-20x20 [--sim] [--compile]", required=False, action="append", help="Plots E and M for a 20x20 lattice")
plots.add_argument(fakeStr("-fakeplotAccFlip"), metavar = "plot-accepted-flips [--sim] [--compile]", required=False, action="append", help="Plots the number of accepted flips as a function of MCCs for T = 1.00, 1.35, 1.70, 2.05 and 2.40")
plots.add_argument(fakeStr("-fakeplotEhist"), metavar = "plot-E-hist [--sim] [--compile]", required=False, action="append", help="Plot energy histogram for temperatures T=1 and T=2.4")


examples = parser.add_argument_group("#Examples#")
examples.add_argument(fakeStr("-ex1"), metavar = "2x2 -T1 8 -filename foo.dat", required=False, action="append", help="Simulates 2x2 lattice from T0=1 to T1=8 and outputs data to foo.dat ")
examples.add_argument(fakeStr("-ex2"), metavar = "20x20 -MCCs 100", required=False, action="append", help="Simulates 20x20 lattice with default values except for MCCs = 100")
examples.add_argument(fakeStr("-ex3"), metavar = "plot-EM-20x20 --sim", required=False, action="append", help="Simulates for default parameters and plots E and M for 20x20 lattice")



def assignDefaults(args, defaults):
  for (arg,val),(darg, dval) in zip(args.__dict__.items(),defaults.items()):
        if (val == None) and (dval != None):
          args.__dict__[arg] = dval
  newargs = dict()
  for key, val in args.__dict__.items():
    if val != None:
      newargs[key] = val
  return newargs

def compile(name):
  os.chdir("../cpp/")
  os.system(f"make {name}")
  os.chdir("../python/")

if __name__ == "__main__":
    if len(sys.argv) == 2:
	    if sys.argv[1] == "compile":
	    	compile("all")
	    	exit()

    args = parser.parse_args()
    defaults = {}
    for arg in args.__dict__.keys():
      defaults[arg] = None

    if args.type=="2x2":
      defaults["T0"] = 1
      defaults["T1"] = 4
      defaults["dT"] = 0.1
      defaults["MCCs"] = int(1e6)
      defaults["filename"] = "../data/2x2EXP.dat"

      args = assignDefaults(args,defaults)

      if args["compile"]:
        compile("2x2")

      subprocess.run(f'../cpp/2x2_lattice.out {args["T0"]} {args["T1"]} {args["dT"]} {args["MCCs"]} {args["filename"]}'.split())

    elif args.type=="20x20":
      defaults["T"] = 1
      defaults["MCCs"] = 10000
      defaults["initSpin"] = 0
      defaults["filename"] = "../data/20x20.dat"
      args = assignDefaults(args,defaults)

      if args["compile"]:
        compile("20x20")
        

      subprocess.run(f'../cpp/20x20_lattice.out {args["T"]} {args["MCCs"]} {args["initSpin"]} {args["filename"]}'.split())

    elif args.type=="parallel":
      defaults["L"] = 20
      defaults["MCCs"] = 1000000
      defaults["saveAfter"] = int(10**4.2)
      defaults["T0"] = 2
      defaults["T1"] = 2.3
      defaults["dT"] = 0.1
      defaults["threads"] = 4
      defaults["filename"] = "parallell_run.dat"

      args = assignDefaults(args, defaults)
      if args["compile"]:
        compile("parallel")

      subprocess.run(f'../cpp/parallel.out {args["L"]} {args["MCCs"]} {args["saveAfter"]} {args["T0"]} {args["T1"]} {args["dT"]} {args["threads"]} {args["filename"]}'.split())

    elif args.type == "plot-exp-2x2":
      defaults["sim"] = False
      defaults["compile"] = False
      args = assignDefaults(args,defaults)
      plotExpValues2x2.main(args["sim"], args["compile"])

    elif args.type == "plot-EM-20x20":
      defaults["sim"] = False
      defaults["compile"] = False
      args = assignDefaults(args,defaults)
      plot20x20.main(args["sim"], args["compile"])
 
    elif args.type == "plot-MCCs-2x2":
      defaults["sim"] = False
      defaults["compile"] = False
      args = assignDefaults(args, defaults)
      plotMCS2x2.main(args["sim"], args["compile"])
    
    elif args.type == "plot-accepted-flips":
      defaults["sim"] = False
      defaults["compile"] = False
      args = assignDefaults(args, defaults)
      plotAcceptedFlips.main(args["sim"], args["compile"])
    
    elif args.type == "plot-E-hist":
      defaults["sim"] = False
      defaults["compile"] = False
      args = assignDefaults(args, defaults)
      plotEnergyHistogram.main(args["sim"], args["compile"])
