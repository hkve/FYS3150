import numpy as np
import matplotlib.pyplot as plt 
import os

sim = True

omegas = [0.01, 0.5, 1, 5]
wf = "T1"
EL = "E0"
filenames = [f"expValues{omega}_{wf}_{EL}.dat" for omega in omegas]
MCCs = 8

E = []
varE = []
r12 = []

if sim:
	for omega, filename in zip(omegas, filenames):
		os.system(f"../compiled/runSingle.exe {wf} {EL} {MCCs} {omega} {filename}")

		E_, EE, r12_ = np.loadtxt(f"../data/{filename}", usecols=(3,4,5)) 

		E.append(E_)
		varE.append(EE-E_*E_)
		r12.append(r12_)

print("omega, E, varE, r12")	
for i in range(len(omegas)):
	print(f"{omegas[i]}, {E[i]}, {varE[i]}, {r12[i]}")
