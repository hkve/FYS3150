import numpy as np
import matplotlib.pyplot as plt

noStable = np.loadtxt("noStable.dat")
Stable = np.loadtxt("Stable.dat")

fig, ax = plt.subplots(nrows=1, ncols=2)
ax[0].set(xscale="log")
ax[1].set(xscale="log")

ax[0].plot(noStable[:,0], noStable[:,1])
ax[1].plot(Stable[:,0], Stable[:,1])

plt.show()