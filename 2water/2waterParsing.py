"""
same as 2H2_H2_dummy_output_parsing_4.py
"""


import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import sys
import get_coords_energy as gce

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
    
def main(fileName, dimension, numbatoms):
    xaxis, VinterSumList, coloumbicEnergySumList, lennardJonesSumList = gce.main(fileName, numbatoms)

    
    VinterSumList = VinterSumList - np.min(VinterSumList)
    print(VinterSumList)

    coloumbicEnergySumList = coloumbicEnergySumList - np.min(coloumbicEnergySumList)
    lennardJonesSumList = lennardJonesSumList - np.min(lennardJonesSumList)


    if dimension == '2': 
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

    elif dimension == 1:
        var1in = input("What is the first variable (x)? ")
        var1 = getVariableIndex(fileName, var1in)
        
        x = arrayVariable(fileName, var1)
        y = PrintArraySCF(fileName)

        editedx = x
        editedy = y - np.min(y)


        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        ax1.scatter(editedx, editedy)
        # ax1.scatter(editedx, VinterSumList)
        ax1.set_xlabel('Oxygen-Oxygen Distance (Angstrom)', fontsize=20)
        ax1.set_ylabel('SCF Energy (Hartree a.u.)', fontsize=20)
        ax1.set_title('H2O H2O Two-Body Energies Inside Clathrate Hydrates', fontsize=25)

        ax2 = fig.add_subplot(212)
        ax2.scatter(editedx, VinterSumList)
        ax2.set_xlabel('Oxygen-Oxygen Distance (Angstrom)', fontsize=20)
        ax2.set_ylabel('SCF Energy (Hartree a.u.)', fontsize=20)
        # ax2.set_title(F'{fileName}')



        # fig = plt.figure()
        # ax1 = fig.add_subplot(411)
        # ax1.scatter(editedx, editedy)
        # # ax1.scatter(editedx, VinterSumList)
        # ax1.set_xlabel(F'{var1in} (Angstrom)')
        # ax1.set_ylabel('SCF Energy (Hartree a.u.)')
        # ax1.set_title('H2O H2O Two-Body Energies Inside Clathrate Hydrates')

        # ax2 = fig.add_subplot(412)
        # ax2.scatter(editedx, VinterSumList)
        # ax2.set_xlabel(F'{var1in} (Angstrom)')
        # ax2.set_ylabel('SCF Energy (Hartree a.u.)')
        # # ax2.set_title(F'{fileName}')

        # ax3 = fig.add_subplot(413)
        # ax3.scatter(editedx, coloumbicEnergySumList)
        # ax3.set_xlabel(F'{var1in} (Angstrom)')
        # ax3.set_ylabel('SCF Energy (Hartree a.u.)')
        # # ax2.set_title(F'{fileName}')

        # ax4 = fig.add_subplot(414)
        # ax4.scatter(editedx, lennardJonesSumList)
        # ax4.set_xlabel(F'{var1in} (Angstrom)')
        # ax4.set_ylabel('SCF Energy (Hartree a.u.)')
        # # ax2.set_title(F'{fileName}')
        
        plt.show()

if __name__ == '__main__':
    print("The arguments are: " + str(sys.argv))


    if len(sys.argv) < 2:
        print('Uh Oh! You did not provide a command line argument that had the name of the file!')
        exit()
        
    else:
        fileName = sys.argv[1]
        print("reading " + fileName)
        
    
    outputvariables, outputvariableslength = variables(fileName)
    print(F'Given these variables: {outputvariables}')
    dimension = input("How many variables do you want to analyze? ")
    dimension = int(dimension)
    numbatoms = input("How many atoms are in your input besides the dummy atom? ")
    numbatoms = int(numbatoms)

    main(fileName, dimension, numbatoms)