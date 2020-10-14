import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys

filename = sys.argv[1]
dt = float(sys.argv[2])
N = int(sys.argv[3])
method = int(sys.argv[4])
beta = float(sys.argv[5])
GR = int(sys.argv[6])

# filename = "sys1"
# dt = 0.1
# N = 10000
# method = 1
# beta = 2
# GR = 1

subprocess.run(f'./main.exe {filename} {dt} {N} {method} {beta} {GR}'.split())
with open(f"data/{filename}.dat", "r+") as file:
    data = file.readlines()
planets = 4
dat0 = []
for i in range(planets):
    dat0.append(data[i*8+1])
    dat0.append(data[i*8+2])

dat = []
for c in dat0:
    dat.append(np.array( list([float(a) for a in c.rstrip("\n").split(",")[:-1]]) ))
c = dat
for i in range(planets):
    plt.plot(np.array(c[i*2]) , np.array(c[i*2+1]), label=str(i) )

plt.grid()
plt.legend()
plt.title(f"dt: {dt}, N: {N}, method: {method}, beta: {beta}, GR: {GR}")
plt.show()