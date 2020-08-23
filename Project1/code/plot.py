import numpy as np
import matplotlib.pyplot as plt

filename = "sol_data"


for i in range(1,4):
	read_filename = filename + str(i) + ".txt"
	data = np.loadtxt(read_filename, delimiter=",")
	plt.plot(data[:,0], data[:,1], label="n = 10e"+str(i))
	x = data[:,0]

analytical = 1-(1-np.exp(-10))*x - np.exp(-10*x)

plt.plot(x, analytical, label="analytical")
plt.legend()
plt.show()