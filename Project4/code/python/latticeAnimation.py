import matplotlib.pyplot as plt
import numpy as np

# Super quick, bytt L til å matche det som er i c++
data = np.loadtxt("../data/lattice.out")

fix, ax = plt.subplots()
L = 300 

while True:
	for i in range(data.shape[0]):
		data_ = data[i].reshape(L,L)
		ax.cla()
		ax.imshow(data_, interpolation="nearest")
		plt.pause(0.5)
plt.show()

