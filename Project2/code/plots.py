import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os as os
import sys
import platform

def get_flags():
	flags = []
	try:
		if sys.argv[1][0] == "-": # check to see whether any flag is given
			for flag in sys.argv[:][1:]:
				flags.append(flag[:][1:])
		else:
			raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")
	except IndexError:
		raise Exception("A flag was not given, nothing was executed. Try the flag -h for a list of available flags and operations.")

	return flags

def check_compile():
	"""
	Function for compiling needed programs in case of missing .exe files
	"""
	files = os.listdir() # Gets filenames from directory 

	if not "data" in files: # Check for missing directory to store datafiles
		os.mkdir("data")
		print("Created missing data directory")

	program_names = ["main.exe", "BucklingBeam.exe"] # Program names
	program_makefile_names = ["compile_main", "compile_bb"] # Program names in the makefile

	for i, name in enumerate(program_names): # Loops over all programs that should be compiled
		if not name in files: # If missing, compile it
			os.system("make " + program_makefile_names[i]) 
			print(f"Compiled {program_names[i]}")

def prase_flags(flags):
	if "h" in flags:
		print("The avalaible flags are")
		print("Not implemented yet...")

if __name__ == "__main__":
	
	flags = get_flags()
	
	check_compile()

	prase_flags(flags)