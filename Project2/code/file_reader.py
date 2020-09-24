import numpy as np

class Run():
	"""
	Class for holding the data for each run. Dosen't do anything at the moment but migth be usefull? 
	at least we don't need to deal with 3D arrays. 
	"""
	def __init__(self, n_iter, N, vals, vecs):
		self.n_iter = n_iter
		self.N = N
		self.vals = vals
		self.vecs = vecs

	def sort_vals_and_vecs(self):
		sort = np.argsort(self.vals) # Indexes of sort after eigenvalue

		self.vals = self.vals[sort] # Swap values
		self.vecs = self.vecs[:,sort] # Swap vectors



def read_data_file(filename):
	"""
	Definitely one of the worst functions I have written, but it works. 
	Im sorry Carl

	args: Take a filename for one of the problems
	returns: A list of instances of Run, such that they can be easily handled 

	Good luck
	"""
	prop_line = True
	runs = []
	with open(filename) as file:
		for line in file:
			if prop_line == True:
				line_counter = 0
				props = line.split()
				n_iter = int(props[0])
				N = int(props[1])

				prop_line = False

				vals = np.zeros(N, dtype=float)
				vecs = np.zeros((N,N), dtype=float)

			else:
				line = line.split()
				
				for j in range(len(line)):
					if j == 0:
						vals[line_counter] = float(line[0])
					else:
						vecs[line_counter,j-1] = float(line[j])

				if line_counter == N-1:
					prop_line = True
					run = Run(n_iter, N, vals, vecs)
					run.sort_vals_and_vecs()
					runs.append(run)
					continue

				line_counter += 1

	return runs 

if __name__ == "__main__":
	# Usage, call the function and enter a filename (no data/)
	runs = read_data_file("data/BucklingBeam.dat")
	
	# Now we can iterate over each run and get som props, ie
	for run in runs:
		print(run.n_iter, run.N)
		print(run.vals)
		print(run.vecs)
	
