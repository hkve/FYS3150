import numpy as np
import matplotlib.pyplot as plt

# Just hardcoded since python will eventually run C++, this must match what you eneter in variateAlphaNoCoulomb.cpp
alphaStart = 0.5
alphaEnd = 1.5
dAlpha = 0.1

N_alpha = round(((alphaEnd-alphaStart)/dAlpha),0)
N_alpha = int(N_alpha) + 1


alphas = np.linspace(alphaStart, alphaEnd, N_alpha, endpoint=True)

fig, axes = plt.subplots(nrows=1, ncols=2)
axes[0].set(xscale="log")
axes[1].set(xscale="log")
for i, alpha in enumerate(alphas):
	data = np.loadtxt(f"../data/varAlphaNoCoulomb{i}.dat")
	cyles, E, E2 = data[:,0], data[:,1], data[:,2]

	varE = E2-E*E
	skip = 0
	axes[0].plot(cyles[skip:], E[skip:], label=fr"$\alpha$ = {alpha:.1f}")
	axes[1].plot(cyles[skip:], varE[skip:], label=fr"$\alpha$ = {alpha:.1f}")

axes[0].legend()
axes[1].legend()
plt.show()