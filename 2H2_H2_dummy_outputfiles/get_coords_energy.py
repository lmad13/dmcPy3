"""
Jordyn Kim

this files pulls out the coordinates and the energy for each energy calculation for scan files and stores it in a xyz file
"""


import sys

fileName = sys.argv[1]
fileName2 = sys.argv[2]
numbatoms = input("How many atoms are in your input besides the dummy atom? ")
numbatoms = int(numbatoms)


def get_SCF(outputFile):
        readFile = open(outputFile,'r')

        SCFlist = []

        for line in readFile:
                if 'SCF Done' in line:
                        line = line.split()
                        SCFlist.append(line[4])

        readFile.close()
        return SCFlist


def get_coords(outputFile, writeFile, SCFlist):
        readFile = open(outputFile,'r')
        writeFile = open(writeFile, 'w+')
        SCFlist = SCFlist

        elementLookUp = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'X']
        countdown = 0
        scfcount = 0 
        emptywrite = ''
        

        for line in readFile:
                if "Standard orientation:" in line:
                        countdown = 6 + numbatoms
                        print("New Run")
                        writeFile.write('{} \n'.format(numbatoms))
                        writeFile.write(' {} \n'.format(SCFlist[scfcount] ))
                        scfcount += 1
                        atomcount = 0
                if countdown > numbatoms:
                        countdown = countdown - 1
                if countdown <= numbatoms and countdown > 0:
                        linelist = line.split()
                        elementnumb = int(linelist[1])
                        element = elementLookUp[elementnumb]
                        linelist[1] = element
                        print(linelist)
                        writetxt = ' {} \t {} \t {} \t {} \n'.format(linelist[1], linelist[3], linelist[4], linelist[5])
                        writeFile.write(writetxt)
                        countdown = countdown - 1

        readFile.close()
        writeFile.close()
        

SCFlist = get_SCF(fileName)
get_coords(fileName, fileName2, SCFlist)


			
