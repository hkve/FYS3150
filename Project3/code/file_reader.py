from getInitialConditions import getUUIDs

class Body():
	def __init__(self, m):
		self.m = m

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
		dt = float
		
		name = get_key(UUID, UUIDs)
		m = Masses[name]
		system[name] = Body(m)

		
		for i, line in enumerate(file):
			if info_line == True:
				print(line)
				"""
				UUID = int(line.split(",")[0])
				name = get_key(UUID, UUIDs)
				print(name)
				"""
			if (i+1)%new_body == 0:
				info_line = True
			else: 
				info_line = False
		
read_data_file("dump.dat")
