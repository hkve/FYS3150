from getInitialConditions import getUUIDs
import numpy as np

class Body():
	def __init__(self, m, UUID, N):
		self.m = m
		self.UUID = UUID
		self.r = np.zeros((3, N))
		self.v = np.zeros((3, N))

def get_key(UUID, UUIDs): 
    for key, value in UUIDs.items(): 
         if UUID == value: 
             return key 

def read_data_file(filename):
	system = {}
	SAVED_DIR = "data/"

	info_line = True

	skipper = 0
	with open(SAVED_DIR + filename) as file:
		UUIDs, Masses = getUUIDs("initData/bodyUUID.txt")
		
		skipper = 1
		for i, line in enumerate(file):
			LINE = line.rstrip().split(",")
			
			if info_line == True:
				UUID = int(LINE[0])
				name = get_key(int(LINE[0]), UUIDs)
				m = Masses[name]

				system["dt"] = float(LINE[1])
				system["N"] = int(LINE[2])
				system["method"] = int(LINE[3])
				system["N_write"] = int(LINE[4])

				system[name] = Body(m, UUID, system["N_write"])

				info_line = False
			
			elif LINE[0] == "*":
				info_line = True
				skipper += 8
			else:
				k = i-skipper			
				for j in range(len(LINE[:-1])):
					if k < 3:
						system[name].r[k, j] = float(LINE[j])
					else:
						system[name].v[k-3, j] = float(LINE[j])

	return system
		
if __name__ == "__main__":
	"""
	Enter filename of data file (stored in data folder) and a dictionary containg all information is returned
	"""
	# If file has has Sun, Earth, Jupiter
	system = read_data_file("SunEarthStable.dat")
	print(len(system["Earth"].r[0]))
	"""
	"""
