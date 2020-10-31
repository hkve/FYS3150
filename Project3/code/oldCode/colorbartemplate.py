import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from colour import Color

#eksampel/template med gauss-dist

NoOfColors = 10
colors = list(Color("cyan").range_to(Color("orange"),NoOfColors)) 
colors = [color.get_rgb() for color in colors]


#dette er verdiene som blir langs colorbaren.
sigs = np.linspace(1,10, NoOfColors, endpoint=True)

f, ax = plt.subplots(1,1)
x = np.linspace(-10, 10, 1000)
def gauss(sig, x):
    return 1/np.sqrt(2*np.pi*sig**2)*np.exp(-1/2*x**2/sig**2)
for i,sig in enumerate(sigs):
    ax.plot(x, gauss(sig, x), color = colors[i], alpha=0.7)

cmap = mpl.colors.ListedColormap(colors)

norm = mpl.colors.Normalize(vmin = sigs[0],vmax = sigs[-1]) # vmin er den verdien som tilsvarer farge #1 (cyan i dette tilfellet) og vmax er den verdien som tilsvarer den siste fargen (oransje idette tilfallet)

cbar = ax.figure.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, orientation='vertical',ticks=sigs )
cbar.ax.set_ylabel(r'   $\sigma$', rotation=0)
plt.show()
