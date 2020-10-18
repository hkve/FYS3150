import numpy as np
from file_reader import read_data_file
from getInitialConditions import setInitialConditions as sic
import matplotlib.pyplot as plt
from subprocess import run


G_ = 6.67408e-11
AU = 1.496e+11
ME = 5.972e+24
G = G_/AU**3*(3600*24*365)**2*ME
c = 299792458/AU*3600*24*365

body_dict = {"Sun": [0,0,0,0,0,0], "Mercury": [0.3075,0,0,0,12.44,0]}
dt = 0.00001
N = 10000000
dpts =  N

sic("GR.dat", body_dict)

run(f'python3.8 master.py {dt} {N} -sys initData/GR.dat -out GR2.out -dpts {dpts} -time years --GR'.split() )
run(f'python3.8 master.py {dt} {N} -sys initData/GR.dat -out notGR2.out -dpts {dpts} -time years'.split() )

systemGR = read_data_file(f"GR2.out")
system = read_data_file(f"notGR2.out")

r = system["Mercury"].r.T
rs = system["Sun"].r.T # r sol
rrel = r-rs #relativ

rGR = systemGR["Mercury"].r.T
rsGR = system["Sun"].r.T
rrelGR = rGR-rsGR


dist = np.linalg.norm(rrel, axis=1)
distGR = np.linalg.norm(rrelGR,axis=1)
print(dist)
print()
print(distGR)

idx = [] # skal fylles med indeks til periheler
idxGR = [] # --""--
#Dette er min metode for å finne minimumspunktene. den funker sånn ish. dersom den gir feil verdier, så kan en finne indeks manuelt ved å se på plottet av dist og manuelt sette nideksen til siste perihel
ismin = np.r_[True, dist[1:] < dist[:-1]] & np.r_[dist[:-1] < dist[1:], True] #noe jeg fant på stackoverflow. 
isminGR = np.r_[True, distGR[1:] < distGR[:-1]] & np.r_[distGR[:-1] < distGR[1:], True]
for i, a in enumerate(ismin):
    if a:
        idx.append(i)
for i, a in enumerate(isminGR):
    if a:
        idxGR.append(i)

print(len(idx), len(idxGR))
theta_p = np.rad2deg(np.arctan(rrel[idx[-1],1]/rrel[idx[-1],0]))
theta_pGR = np.rad2deg(np.arctan(rrelGR[idxGR[-1],1]/rrelGR[idxGR[-1],0]))

print(theta_p*3600,theta_pGR*3600) # arcseconds i perihel
print(idx[-1], idxGR[-1])
plt.plot(distGR, label="distGR")
plt.plot(dist, label="dist")
plt.legend()
plt.show()
plt.plot(rGR[:,0], rGR[:,1], label="distGR")
plt.plot(r[:,0], r[:,1], label="dist")
plt.legend()
plt.show()
  
    
