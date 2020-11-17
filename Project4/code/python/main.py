import numpy as np
import sys, argparse, subprocess, argcomplete


import plot20x20, plotExpValues2x2, plotMCS2x2, plotAcceptedFlips


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

parser.add_argument('type',type=str,  choices = ["2x2", "20x20", "plot-exp-2x2", "plot-EM-20x20", "plot-MCCs-2x2", "plot-accepted-flips"],help="Follow program name by optional program specific parameters listed below")
parser.add_argument("-MCCs", metavar="",type=int, help="int. Number of MCCs to perform")
parser.add_argument("-T", metavar="",type=float, help="float. System temperature")
parser.add_argument("-T0", metavar = "", type=float, help="float. Start temperature")
parser.add_argument("-T1", metavar = "", type=float, help="float. End temperature")
parser.add_argument("-dT", metavar = "", type=float, help="float. Temperature step")
parser.add_argument("-initSpin", metavar="", type=int, help="int. Init spin configuration. 1: all up, -1: all down, 0: random")
parser.add_argument("-filename", metavar="", type=str, help="str. Output filename")
parser.add_argument("--sim", action="store_true", help="Re-run simulations in plotting code with default parameters. May take time")


progparams = parser.add_argument_group("Program parameters (specifics)")
progparams.add_argument(fakeStr("-fake2x2"), metavar = "2x2 [-T0] [-T1] [-dT] [-MCCs] [-filename]", required=False, action="append", help="2x2 lattice from temperature T0 to T1 with step dT. Defaults: T0=1, T1=4, dT=0.1, MCCs=100000, filename=../data/2x2EXP.out")
progparams.add_argument(fakeStr("-fake20x20"), metavar = "20x20 [-T] [-MCCs] [-initSpin] [-filename]", required=False, action="append", help="20x20 lattice at temp. T. Defaults: T=1, MCCs=100000, initSpin = 0, filename=../data/20x20.out")
progparams.add_argument(fakeStr("-fakeplot2x2"), metavar = "plot-exp-2x2 [--sim]", required=False, action="append", help="Plots exp. values for a 2x2 lattice")
progparams.add_argument(fakeStr("-fakeplot2x2_MCCs"), metavar = "plot-MCCs-2x2 [--sim]", required=False, action="append", help="Plots exp. values vs MCCs for a 2x2 lattice")
progparams.add_argument(fakeStr("-fakeplot20x20"), metavar = "plot-EM-20x20 [--sim]", required=False, action="append", help="Plots E and M for a 20x20 lattice")
progparams.add_argument(fakeStr("-fakeplotAccFlip"), metavar = "plot-accepted-flips [--sim]", required=False, action="append", help="DESCRIPTION MISSING")


examples = parser.add_argument_group("#Examples#")
examples.add_argument(fakeStr("-ex1"), metavar = "2x2 -T1 8 -filename ../data/foo.dat", required=False, action="append", help="Simulates 2x2 lattice from T0=1 to T1=8 and outputs data to foo.dat ")
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


if __name__ == "__main__":
    args = parser.parse_args()
    defaults = {}
    for arg in args.__dict__.keys():
      defaults[arg] = None
   

    if args.type=="2x2":
      defaults["T0"] = 1
      defaults["T1"] = 4
      defaults["dT"] = 0.1
      defaults["MCCs"] = int(1e6)
      defaults["filename"] = "../data/2x2EXP.out"

      args = assignDefaults(args,defaults)
      subprocess.run(f'../cpp/2x2_lattice.out {args["T0"]} {args["T1"]} {args["dT"]} {args["MCCs"]} {args["filename"]}'.split())

    elif args.type=="20x20":
      defaults["T"] = 1
      defaults["MCCs"] = 10000
      defaults["initSpin"] = 0
      defaults["filename"] = "../data/20x20.out"
      args = assignDefaults(args,defaults)

      subprocess.run(f'../cpp/20x20_lattice.out {args["T"]} {args["MCCs"]} {args["initSpin"]} {args["filename"]}'.split())

    elif args.type == "plot-exp-2x2":
      defaults["sim"] = False
      args = assignDefaults(args,defaults)
      plotExpValues2x2.main(args["sim"])

    elif args.type == "plot-EM-20x20":
      defaults["sim"] = False
      args = assignDefaults(args,defaults)
      plot20x20.main(args["sim"])
 
    elif args.type == "plot-MCCs-2x2":
      defaults["sim"] = False
      args = assignDefaults(args, defaults)
      plotMCS2x2.main(args["sim"])
    
    elif args.type == "plot-accepted-flips":
      defaults["sim"] = False
      args = assignDefaults(args, defaults)
      plotAcceptedFlips.main(args["sim"])