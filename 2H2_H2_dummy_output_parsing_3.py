"""
Since the fidget scanner varies both A6 and A7, I only get the data from scans where the two H connected to the dummy atom are in line
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

def PrintArray(fileName):
    readFile = open(fileName, 'r')
    data = readFile
    b6 = []
    a6 = []
    a7 = []
    scf = []
    printThisLine = False
    for line in data:
        if 'Summary of the potential surface scan:' in line:
            printThisLine = True
        if 'JKIM22' in line:
            printThisLine = False
        if printThisLine:
            listline = line.split()
            if len(listline) == 5:
                if len(listline[0]) < 6 & len(listline[1]) > 2:
                    checka6 = 90 - float(listline[2])
                    checka7 = float(listline[3]) - 90
                    if checka6 == checka7:
                        b6.append(float(listline[1]))
                        a6.append(float(listline[2]))
                        a7.append(float(listline[3]))
                        scf.append(float(listline[4]))
    
    readFile.close()
    return np.array(b6), np.array(a6), np.array(a7), np.array(scf)

b6, a6, a7, scf = PrintArray(fileName)

x = b6
y = a6
z = scf

editedx = x - x[0]
editedy = y - y[0] 

fig = plt.figure()
ax1 = fig.add_subplot(111, projection = '3d')
ax1.scatter(editedx, editedy, z, c=z, cmap='RdBu')
ax1.set_xlabel('B6 (Angstrom)')
ax1.set_ylabel('A6 (Degrees)')
ax1.set_zlabel('Energy (a.u.)')
ax1.set_title(F'{fileName}')
plt.show()
