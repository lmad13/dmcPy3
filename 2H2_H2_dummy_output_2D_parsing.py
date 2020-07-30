"""
graph SCF energies for 2H2_H2_dummy_input_0.log and 2H2_H2_dummy_input_1.log
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import sys

print("The arguments are: " + str(sys.argv))
if len(sys.argv) < 2:
    print('Uh Oh! You did not provide a command line argument that had the name of the file!')
    exit()
else:
    fileName = sys.argv[1]
    print("reading " + fileName)

def PrintArraySCF(fileName):
    readFile = open(fileName, 'r')
    scflist = []
    data = readFile
    for line in data:
        if "SCF Done" in line:
            energy_line = line
            words = energy_line.split()
            energy = words[4]
            energy = float(energy)
            scflist.append(energy)
    
    readFile.close()
    return np.array(scflist)

def PrintArrayb6(fileName):
    readFile = open(fileName, 'r')
    data = readFile
    variablelist = []
    printThisLine = False
    for line in data:
        if 'Summary of the potential surface scan:' in line:
            printThisLine = True
        if 'JKIM22' in line:
            printThisLine = False
        if printThisLine:
            listline = line.split()
            if len(listline) == 3:
                if len(listline[0]) < 6 & len(listline[1]) > 2:
                    b6 = float(listline[1])
                    variablelist.append(b6)
    
    readFile.close()
    return np.array(variablelist)

x = PrintArrayb6(fileName)
y = PrintArraySCF(fileName)

editedx = x - x[0]
editedy = y - np.min(y)

editedxright = -1*editedx

combinex = np.concatenate((editedx, editedxright[1:]))
combiney = np.concatenate((editedy, editedy[1:]))

print("B6: ", combinex)
print("SCF energy: ", combiney)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(combinex,combiney)
ax1.set_xlabel('B6 (Angstrom)')
ax1.set_ylabel('Energy (a.u.)')
ax1.set_title(F'{fileName}')
plt.show()