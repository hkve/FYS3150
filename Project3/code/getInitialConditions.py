from astroquery.jplhorizons import Horizons
import os as os

UIDs = { # Fant ikke en bedre måte å gjøre dette på
	"Sun": 10,
	"Mercury": 199,                                                              
	"Venus": 299,                                                             
	"Moon": 301,                                                             
	"Earth": 399,                                                       
	"Mars": 499,                                                                 
	"Phobos": 401,                                                             
	"Deimos": 402,                                                           
	"Io": 501,                                                                  
	"Europa": 502,                                                            
	"Ganymede": 503,                                                         
	"Callisto": 504,
	"Jupiter": 599,
	"Saturn": 699,  
	"Uranus": 799, 
	"Neptune": 899,  
	"Pluto": 999     
  	}

Masses = {
	"Sun": 1.989E30,
	"Mercury": 3.285E23,                                                              
	"Venus": 4.867E24,                                                             
	"Moon": 7.34767309E22,                                                             
	"Earth": 5.972E24,                                                       
	"Mars": 6.39E23,                                                                 
	"Phobos": 1.08e16,                                                             
	"Deimos": 2.0e15,                                                           
	"Io": 8.9319e22,                                                                  
	"Europa": 4.799844e22,                                                            
	"Ganymede": 1.4819e23,                                                         
	"Callisto": 1.075938e23,
	"Jupiter": 1.898e27,
	"Saturn": 5.683e26,  
	"Uranus": 8.681e25, 
	"Neptune": 1.024e26,  
	"Pluto": 1.30900e22  
}  	

def grabBody(bodyName, date=None):
	"""
	Takes a body name from UIDs dictionary and returns initial position and velocity
	
	Args: 
		bodyName: String, name of the desired body, example "Sun"
		date: 	  String, YYYY-MM-DD, Date for intial conditions, default is todays date

	Returns:
		init:     String ready to be written to file
	"""
	UID = UIDs[bodyName]
	obj = Horizons(id=UID, location="500@0", epochs=date,  id_type="majorbody")
	v = obj.vectors()
	
	init = [str(v[elm][0]) for elm in ["x", "y", "z", "vx", "vy", "vz"]]
	m = Masses[bodyName]
	init.append(str(m))
	init.insert(0, str(UID))
	init = ",".join(init) 
	
	return init + "\n"

def getInitialCondition(filename, bodies, date=None):
	"""
	Takes a list of bodies and finds initial conditions

	Args:
		filename: String, name for file to write to
		bodies:   String list, name of bodies for initial conditions
		date: 	  String, YYYY-MM-DD, Date for intial conditions, default is todays date

	"""
	DIR = os.listdir()
	if not "data" in DIR: # Check for missing directory to store datafiles
		os.mkdir("data")
	
	filename = "data/" + filename

	with open(filename, "w+") as file:
		for body in bodies:
			if body not in UIDs.keys():
				print(f"{body} not available, skips")
				continue

			file.write(grabBody(body, date))

if __name__ == "__main__":
	# Examples
	filename = "SunEarthJupiter_init.dat"
	bodies = ["Sun", "Earth", "Jupiter", "Moren din"]
	getInitialCondition(filename, bodies)
	
	# All planets, Earth moon, Mars moons and Jupiters 4 most massive moons (and Pluto <3)
	filename = "SolarSystem_init.dat"
	getInitialCondition(filename, UIDs.keys())