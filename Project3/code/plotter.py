import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys
dt = float(sys.argv[1])
N = int(sys.argv[2])
method = int(sys.argv[3])
subprocess.run(f"./main.out sys1.txt {dt} {N} {method}".split())
with open("data/TEMPNAME.dat", "r+") as file:
    data = file.readlines()
planets = 3
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
plt.title(f"dt: {dt}, N: {N}, method: {method}")
plt.show()