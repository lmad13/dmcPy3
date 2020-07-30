"""
combine 2H2_H2_dummy_input_5.log and 2H2_H2_dummy_input_4.log to predict what would happen to the PES if the dummy atom from 2H2_H2_dummy_input_5.log approached the other water
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt

def PrintArraySCF(fileName):
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
                    scf = float(listline[3])
                    variablelist.append(scf)
    
    readFile.close()
    return np.array(variablelist)
                  

scfarray = PrintArraySCF('2H2_H2_dummy_input_5.log')
scfarray2 = PrintArraySCF('2H2_H2_dummy_input_4.log')

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
                  
b6array = PrintArrayb6('2H2_H2_dummy_input_5.log')
b6array2 = PrintArrayb6('2H2_H2_dummy_input_4.log')

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

a5array = PrintArraya5('2H2_H2_dummy_input_5.log')
a5array2 = PrintArraya5('2H2_H2_dummy_input_4.log')

x = np.concatenate( (b6array, b6array2))
y = np.concatenate( (a5array, a5array2))
z = np.concatenate( (scfarray, scfarray2))

fig = plt.figure()
ax1 = fig.add_subplot(111, projection = '3d')
ax1.scatter(x, y, z)
ax1.set_xlabel('B6 (Angstrom)')
ax1.set_ylabel('A5 (Degrees)')
ax1.set_zlabel('Energy (a.u.)')
plt.show()

print(x)
print(y)
print(z)

print(np.min(z))