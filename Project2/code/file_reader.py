import numpy as np

class Run():
	"""
	Class for holding the data for each run. Dosen't do anything at the moment but migth be useful? 
	at least we don't need to deal with 3D arrays. 
	"""
	valid_properties = ['n_iter', 'N', 'rho_max', 'omega_r'] # The run-dependent properties will be stored as a dict
	def __init__(self, vals, vecs, **kwargs):
		self.properties = {} # Depending on what type of run it is, different properties are stored. (Could potentially make subclasses)
		for key, item in kwargs.items():
			if key in self.valid_properties:
				self.properties[key] = item
			else:
				raise Exception(f"Could not understand argument '{key}'. Valid properties as {valid_properties}")
		self.vals = vals
		self.vecs = vecs

	def sort_vals_and_vecs(self):
		sort = np.argsort(self.vals) # Indexes of sort after eigenvalue

		self.vals = self.vals[sort] # Swap values
		self.vecs = self.vecs[:,sort] # Swap vectors

	def __call__(self, prop):
		# Calling the instance of the class with a string denoting the value wanted returns the property correspondingly.
		if prop in self.properties:
			return self.properties[prop]
		else:
			raise Exception(f"The given property is not recorded in this run. ('{prop}')")




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
	valid_props = ['n_iter', 'N', 'rho_max', 'omega_r'] # valid properties for a run
	with open(filename) as file:
		for line in file:
			if prop_line == True:
				line_counter = 0
				prop_vals = line.split()

				run_properties = {} # for storing relevant properties of the particular run
				for i, value in enumerate(prop_vals):
					run_properties[valid_props[i]] = eval(value)
				N = run_properties['N']

				prop_line = False

				vals = np.zeros(N, dtype=float)
				vecs = np.zeros((N,N), dtype=float)

			else:
				line = line.split()
				
				for j in range(len(line)):
					if j == 0:
						# first value is the eigenvalue
						vals[line_counter] = float(line[0])
					else:
						# the following are the eigenvectors
						vecs[line_counter,j-1] = float(line[j])

				if line_counter == N-1:
					prop_line = True
					run = Run(vals, vecs, **run_properties)
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
	
