"""
same as 2H2_H2_dummy_output_parsing.py, but I will try to get variables from command line input
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
    
    
def variables(fileName):
    readFile = open(fileName, 'r')
    data = readFile
    printThisLine = False
    for line in data:
        if 'Summary of the potential surface scan:' in line:
            printThisLine = True
        if '-----------' in line:
            printThisLine = False
        if printThisLine:
            listline = line.split()
            if 'SCF' in listline:
                variables = listline[1:-1]
                linelength = len(listline)

    readFile.close()
    return variables, linelength

outputvariables, outputvariableslength = variables(fileName)
print(F'Given these variables: {outputvariables}')
dimension = input("How many variables do you want to analyze? ")


if dimension == '2': 
    
    def variables(fileName):
        readFile = open(fileName, 'r')
        data = readFile
        printThisLine = False
        for line in data:
            if 'Summary of the potential surface scan:' in line:
                printThisLine = True
            if '-----------' in line:
                printThisLine = False
            if printThisLine:
                listline = line.split()
                if 'SCF' in listline:
                    variables = listline[1:-1]   

        readFile.close()
        return variables
    
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

    def getVariableIndex(fileName, variable):
        readFile = open(fileName, 'r')
        var = variable
        data = readFile
        printThisLine = False
        for line in data:
            if 'Summary of the potential surface scan:' in line:
                printThisLine = True
            if '-----------' in line:
                printThisLine = False
            if printThisLine:
                listline = line.split()
                if var in listline:
                    index = listline.index(var)    

        readFile.close()
        return index

    def arrayVariable(fileName, variable):
        readFile = open(fileName, 'r')
        var = variable
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
                if len(listline) == outputvariableslength:
                    if len(listline[0]) < 4 & len(listline[1]) > 2:
                        desiredvar = float(listline[var])
                        variablelist.append(desiredvar)

        readFile.close()
        return np.array(variablelist)

    var1in = input("What is the first variable (x)? ")
    var2in = input("What is the second variable (y)? ")

    var1 = getVariableIndex(fileName, var1in)
    var2 = getVariableIndex(fileName, var2in)
    
    x = arrayVariable(fileName, var1)
    y = arrayVariable(fileName, var2)
    z = PrintArraySCF(fileName)

    editedx = x
    editedy = y

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    ax1.scatter(editedx, editedy, z, c=z, cmap='RdBu')
    ax1.set_xlabel(F'{var1in} (Angstrom)')
    ax1.set_ylabel(F'{var2in} (Degrees)')
    ax1.set_zlabel('Energy (a.u.)')
    ax1.set_title(F'{fileName}')
    plt.show()
    


    
