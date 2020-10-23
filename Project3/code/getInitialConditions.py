from astroquery.jplhorizons import Horizons
import os as os

def getUUIDs(filename):
	"""
	Reads and stores names, UUIDs and masses of bodies

	Args:
		filename: String, filename where data is stored

	Returns:
		UUIDs: Dictionary with names as keys and UUID as value
		Masses: Dictionary with names as keys and mass as value 
	"""

	UUIDs = {}
	Masses = {}

	with open(filename) as file:
		for line in file:
			LINE = line.strip().split(",")
			name, UUID, m = LINE
			UUIDs[name] = int(UUID)
			Masses[name] = float(m)

	return UUIDs, Masses

def grabBody(UUID, date=None):
	"""
	Takes a UUID and returns string of UUID, initial position and velocity
	
	Args: 
		UUID: 	  The unique id for a body
		date: 	  String, YYYY-MM-DD, Date for intial conditions, default is todays date

	Returns:
		init:     String of UUID and initial conditions (UUID,x,y,z,vx,vy,vz)
	"""

	obj = Horizons(id=UUID, location="500@0", epochs=date,  id_type="majorbody")
	v = obj.vectors()
	
	init = [str(v[elm][0]) for elm in ["x", "y", "z", "vx", "vy", "vz"]]
	init.insert(0, str(UUID))
	init = ",".join(init) 
	
	return init

def getInitialCondition(filename, bodies=None, date=None):
	"""
	Takes a list of bodies and writes initial condtions to filename

	Args:
		filename: String, name for file to write to
		bodies:   String list, name of bodies for initial conditions. If None all bodies are written
		date: 	  String, YYYY-MM-DD, Date for intial conditions, default is todays date

	"""
	
	DIR2SAVE = "initData"
	UUIDs, Masses  = getUUIDs(DIR2SAVE+"/bodyUUID.txt")

	DIR = os.listdir()

	if not DIR2SAVE in DIR: # Check for missing directory to store datafiles
		os.mkdir(DIR2SAVE)
	
	filename = DIR2SAVE + "/" + filename

	if bodies == None:
		bodies = UUIDs.keys() # If None is given, wrirte alle bodies

	n_bodies = len(bodies)
	with open(filename, "w+") as file:
		for i, body in enumerate(bodies):
			if body not in UUIDs.keys():
				print(f"{body} not available, skips")
				continue

			UUID = UUIDs[body]
			str2write = grabBody(UUID) # Gets x and v values
			str2write =  str2write + "," + str(Masses[body]/Masses["Earth"]) 

			if i < n_bodies-1:
				str2write += "\n"
			
			file.write(str2write)

	print(f"Saved to {filename}, wrote {len(bodies)} bodies")

def setInitialConditions(filename, body_dict, fixedCoM = False):
	"""
	For manuel setting of initial conditions

	Args:
		filename: String, name of file to write to
		body_dict: Dictionary, {bodyname: [x,y,z,vx,vy,vz]}
		fixedCoM: Bool, whether to adjust the initial conditions such that the center of mass has no velocity
	"""

	DIR2SAVE = "initData"
	UUIDs, Masses  = getUUIDs(DIR2SAVE+"/bodyUUID.txt")


	for body in body_dict:
		init_len = len(body_dict[body])
		if init_len != 6: 
			print(f"Body {body} has wrong init condtions, excepted 6 got {init_len}")
			exit()

		for i in range(6):
			body_dict[body][i] = str(body_dict[body][i])


	if fixedCoM:
		# Disclaimer: this is slightly messy, but it works
		P = [0,0,0] # total momentum of system
		R = [0,0,0]
		Mtot = 0 # total mass of system
		for body, values in body_dict.items():
			M = Masses[body] # mass of current body
			P = [p + M * float(v) for p, v in zip(P, values[3:])] # adding momentum of current body to total
			R = [r + M*float(x) for r, x in zip(R, values[0:3])]
			Mtot += M
		V = [p/Mtot for p in P] # center of mass velocity
		R = [r/Mtot for r in R]
		
		for body, values in body_dict.items():
			# adjusting every initial velocity with the CoM velocity
			values[0:3] = [str(float(value) - r) for value, r in zip(values[0:3], R)]
			values[3:] = [str(float(value) - v) for value, v in zip(values[3:], V)]


	dict_len = len(body_dict)
	with open(DIR2SAVE+ "/" + filename, "w+") as file:
		for i, body in enumerate(body_dict):
			line = ",".join(body_dict[body])
			line = str(UUIDs[body]) + "," + line + "," +str(Masses[body]/Masses["Earth"]) 
			
			if i < dict_len-1:
				line += "\n"

			file.write(line)
			

if __name__ == "__main__":
	# Examples
	filename = "SunEarthJupiter_init.dat"
	bodies = ["Sun", "Earth", "Jupiter"]
	getInitialCondition(filename, bodies)
	"""	
	
	# All planets, Earth moon, Mars moons and Jupiters 4 most massive moons (and Pluto <3)
	filename = "SolarSystem_init.dat"
	getInitialCondition(filename)
	# Manual
	body_dict = {"Sun": [0,0,0,0,0,0],
				 "Earth": [1,0,0,6.28318530718,0,0]}
	setInitialConditions("SunEarth_init.dat", body_dict)
	"""