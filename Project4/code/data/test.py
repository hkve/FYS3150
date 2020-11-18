import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("../data/Energy_count24.dat")
tot = np.sum(data)


plt.hist(data, bins=112)
plt.show()