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

scfarray = PrintArraySCF(fileName)

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
            if len(listline) == 4:
                if len(listline[0]) < 4 & len(listline[1]) > 2:
                    b6 = float(listline[1])
                    variablelist.append(b6)
    
    readFile.close()
    return np.array(variablelist)
                  
b6array = PrintArrayb6(fileName)

def PrintArraya5(fileName):
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
            if len(listline) == 4:
                if len(listline[0]) < 4 & len(listline[1]) > 2:
                    a5 = float(listline[2])
                    variablelist.append(a5)
    
    readFile.close()
    return np.array(variablelist)

a5array = PrintArraya5(fileName)

x = b6array
y = a5array
z = scfarray

editedx = x-x[0]
editedy = y-y[0]

editedxright = -1*editedx[1:]
editedyright = -1*editedy[1:]

fig = plt.figure()
ax1 = fig.add_subplot(111, projection = '3d')
ax1.scatter(editedx, editedy, z, c=z, cmap='RdBu')
ax1.set_xlabel('B6 (Angstrom)')
ax1.set_ylabel('A5 (Degrees)')
ax1.set_zlabel('Energy (a.u.)')
ax1.set_title(F'{fileName}')
plt.show()