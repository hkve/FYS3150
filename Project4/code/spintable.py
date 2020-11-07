import numpy as np
import pandas as pd




ddic = {"s0":[], "E0":[], "s1":[], "E1":[], "s2":[], "E2":[], "s3":[], "E3":[]}
up = "\\uparrow"
down = "\\downarrow"
sdic = {1:up, -1:down}
L = 2
types = ["1111", "1110", "1101", "1011", "0111", "0011", "0110", "1100", "1001", "1010", "0101", "0001", "0010", "0100", "1000", "0000"]

def get_E(spins_):
    neighs = np.array([ [[(0,1), (1,0)], [(0,0),(1,1)] ] , [ [(0,0),(1,1)], [(1,0),(0,1)]  ]   ])
    res = 0
    for i in range(L):
        for j in range(L):
            idx = (i,j)
            res += spins_[idx]*(spins_[tuple(neighs[idx][0])] + spins_[tuple(neighs[idx][1])])
    if res != 0:
        return f"{int(-res)}J"
    return f"{0}"

def spin_str(spins):   
    res = " $\\begin{smallmatrix} "
    res += sdic[spins[0,0]] + " & "
    res += sdic[spins[0,1]] + " \\\\ "
    res += sdic[spins[1,0]] + " & "
    res += sdic[spins[1,1]] + " \\end{smallmatrix}$"
    #print(res)
    return res
#H3fgJmPaNAdjyg
def str_to_spin(string):
    res = np.zeros([L,L])
    strdic = {"1":1, "0":-1}
    res[0,0] = strdic[string[0]]
    res[0,1] = strdic[string[1]]
    res[1,0] = strdic[string[2]]
    res[1,1] = strdic[string[3]]
    return res
    

for i,string in enumerate(types):
    spins = str_to_spin(string)
    print(get_E(spins), spins, np.sum(spins))
    ddic[f"E{int(i//4)}"].append(spin_str(spins))
    ddic[f"s{int(i//4)}"].append(f"{get_E(spins)}")
    #print(int(i//4))

pd.set_option('max_colwidth',200)
df = pd.DataFrame(ddic, copy=True, columns=4)
#print(df)
#print(df.values[0])
#print(df.to_latex(escape=False))
        