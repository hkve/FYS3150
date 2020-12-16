import numpy as np
import sys, argparse, subprocess, os

import optimalStepFitter, plotOptimalStep, plotStability, plotEvarEvsAlpha, plot3DVarBeta, plotVarBetaParallel, plotVirial


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
parser._positionals.title = 'Plots'
parser._optionals.title = 'args'
comp = parser.add_argument_group("compile all programs")
comp.add_argument(fakeStr("-fakecompile"), metavar = "compile", required=False, action="append", help="Compiles all C++ programs")

parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS, help='Show this help message and exit.')
argparse._HelpAction(option_strings=['-h', '--help'], dest='help', default='==SUPPRESS==', help='Show this help message and exit.')

parser.add_argument('type',type=str,  choices = ["stability", "optimal-step", "optimal-step-fit", "E-varE", "3D-beta", "beta-parallel", "virial"],help="Follow program name by optional program specific parameters listed below")
parser.add_argument("-MCCs", metavar="",type=int, help="int. Number of Monte Carlo cycles to perform")
parser.add_argument("-alpha", metavar = "", type=float, help="float. α, variational parameter")
parser.add_argument("-beta", metavar="", type=float, help="float. β, variational parameter")

parser.add_argument("--sim", action="store_true", help="Re-run simulations in plotting code with default parameters. May take time")
parser.add_argument("--compile", metavar="", type=str, help="str. compile individual programs")


progparams = parser.add_argument_group("Program parameters (specifics) \nFrom here you can run the C++ programs")
progparams.add_argument(fakeStr("-fakeRunSingle"), metavar = "runSingle [-T0] [-T1] [-dT] [-MCCs] [-filename] [--compile]", required=False, action="append", help="FILL ME")
progparams.add_argument(fakeStr("-fakeAlpha"), metavar = "variateAlpha [-T] [-MCCs] [-initSpin] [-filename] [--compile]", required=False, action="append", help="FILL ME")
progparams.add_argument(fakeStr("-fakeBeta"), metavar="variateBeta [-L] [-MCCs] [-saveAfter] [-T0] [-T1] [-dT] [threads] [-filename] [--compile]", required=False, action="append", help="FILL ME")

plots = parser.add_argument_group("Program plots \nFrom here you can recreate the plots as seen in the report")
plots.add_argument(fakeStr("-fakeplotStability"), metavar = "stability [--sim]", required=False, action="append", help="Plots exp. values for a 2x2 lattice")
plots.add_argument(fakeStr("-fakeplotoptimalstep"), metavar = "opt-step [--sim]", required=False, action="append", help="rrPlots exp. values2x2 lattice ")
plots.add_argument(fakeStr("-fakeplotoptimalstepfit"), metavar = "opt-step-fit", required=False, action="append", help="Plots E and M for a 20x20 lattice")
plots.add_argument(fakeStr("-fakeplotEvarE"), metavar = "E-varE [--sim]", required=False, action="append", help="Plots the number of accepted flips as a function of MCCs for T = 1.00, 1.35, 1.70, 2.05 and 2.40")
plots.add_argument(fakeStr("-fakeplot3Dbeta"), metavar = "3D-beta [--sim]", required=False, action="append", help="Plot energy histogram for temperatures T=1 and T=2.4")
plots.add_argument(fakeStr("-fakeplotbetaparalell"), metavar = "beta-par [--sim]", required=False, action="append", help="Plots the number of accepted flips as a function of MCCs for T = 1.00, 1.35, 1.70, 2.05 and 2.40")
plots.add_argument(fakeStr("-fakeplotvirial"), metavar = "virial[--sim] ", required=False, action="append", help="Plots the number of accepted flips as a function of MCCs for T = 1.00, 1.35, 1.70, 2.05 and 2.40")


examples = parser.add_argument_group("#Examples#")
examples.add_argument(fakeStr("-ex1"), metavar = "3D-beta --sim", required=False, action="append", help="Simulates 3D alpha/beta plot")
#examples.add_argument(fakeStr("-ex2"), metavar = "20x20 -MCCs 100", required=False, action="append", help="Simulates 20x20 lattice with default values except for MCCs = 100")
#examples.add_argument(fakeStr("-ex3"), metavar = "plot-EM-20x20 --sim", required=False, action="append", help="Simulates for default parameters and plots E and M for 20x20 lattice")



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
	os.chdir("../VMC/")
	os.system(f"make {name}")
	os.chdir("../plotting/")

if __name__ == "__main__":
	if not os.path.isdir("../data/"):
		subprocess.run(f"mkdir ../data".split())

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

	elif args.type == "stability":
		defaults["sim"] = False
		args = assignDefaults(args,defaults)
		plotStability.main(args["sim"])

	elif args.type == "optimal-step":
		defaults["sim"] = False
		args = assignDefaults(args,defaults)
		plotOptimalStep.main(args["sim"])
 
	elif args.type == "optimal-step-fit":
		optimalStepFitter.main()

	elif args.type == "E-varE":
		defaults["sim"] = False
		args = assignDefaults(args, defaults)
		plotEvarEvsAlpha.main(args["sim"])
	
	
	elif args.type == "3D-beta":
		defaults["sim"] = False
		args = assignDefaults(args, defaults)
		plot3DVarBeta.main(args["sim"])

	elif args.type == "beta-parallel":
		defaults["sim"] = False
		args = assignDefaults(args, defaults)
		plotVarBetaParallel.main(args["sim"])

	elif args.type == "virial":
		defaults["sim"] = False
		args = assignDefaults(args, defaults)
		plotVirial.main(args["sim"])

	
