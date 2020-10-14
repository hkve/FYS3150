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
	UUIDs, Masses  = getUUIDs("initData/bodyUUID.txt")

	DIR = os.listdir()
	DIR2SAVE = "initData"

	if not DIR2SAVE in DIR: # Check for missing directory to store datafiles
		os.mkdir(DIR2SAVE)
	
	filename = DIR2SAVE + "/" + filename

	if bodies == None:
		bodies = UUIDs.keys() # If None is given, wrirte alle bodies

	with open(filename, "w+") as file:
		for body in bodies:
			if body not in UUIDs.keys():
				print(f"{body} not available, skips")
				continue

			UUID = UUIDs[body]
			str2write = grabBody(UUID) # Gets x and v values
			str2write =  str2write + "," + str(Masses[body]) + "\n" # writes
			file.write(str2write)

	print(f"Saved to {filename}, wrote {len(bodies)} bodies")

if __name__ == "__main__":
	# Examples
	filename = "SunEarthJupiter_init.dat"
	bodies = ["Sun", "Earth", "Jupiter"]
	getInitialCondition(filename, bodies)
	
	# All planets, Earth moon, Mars moons and Jupiters 4 most massive moons (and Pluto <3)
	filename = "SolarSystem_init.dat"
	getInitialCondition(filename)