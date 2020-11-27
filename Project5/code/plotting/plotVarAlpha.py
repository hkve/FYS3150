import numpy as np
import matplotlib.pyplot as plt

alphas = np.linspace(0.6, 1.4, 5, endpoint=True)

fig, axes = plt.subplots(nrows=1, ncols=2)
for i, alpha in enumerate(alphas):
	data = np.loadtxt(f"../data/varAlphaNoCoulomb{i}.dat")
	cyles, E, E2 = data[:,0], data[:,1], data[:,2]

	varE = E2-E*E

	axes[0].plot(cyles, E, label=fr"$\alpha$ = {alpha:.1f}")
	axes[1].plot(cyles, varE, label=fr"$\alpha$ = {alpha:.1f}")

axes[1].legend()
plt.show()