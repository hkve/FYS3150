from getInitialConditions import getUUIDs
import numpy as np

class Body():
	def __init__(self, m):
		self.m = m

	def set_cords(self, r, v):
		self.r = r
		self.v = v	

def get_key(UUID, UUIDs): 
    for key, value in UUIDs.items(): 
         if UUID == value: 
             return key 

def read_data_file(filename):
	system = {}
	SAVED_DIR = "data/"

	info_line = True
	new_body = 8

	with open(SAVED_DIR + filename) as file:
		UUIDs, Masses = getUUIDs("initData/bodyUUID.txt")
		
		first_line = file.readline().rstrip().split(",")
		UUID = int(first_line[0])
		system["dt"] = float(first_line[1])
		system["N"] = int(first_line[2])

		name = get_key(UUID, UUIDs)
		m = Masses[name]
		system[name] = Body(m)

		is_first_line = False

		skipper = 0

		body_data = np.zeros((6, system["N"]))
		for i, line in enumerate(file):
			LINE = line.rstrip().split(",")

			if is_first_line == True:
				print(name)
				system[name].set_cords(body_data[0:3, :], body_data[3:, :]) 
				
				UUID = int(LINE[0])
				name = get_key(UUID, UUIDs)
				m = Masses[name]
				system[name] = Body(m)
				
				skipper += 8
				is_first_line = False

			elif LINE[0] == "*":
				is_first_line = True
			else:
				k = i-skipper
				for j, value in enumerate(LINE[:-1]):
					body_data[k,j] = float(value)	
	print(system["Sun"].r[0,0])
read_data_file("dump.dat")
