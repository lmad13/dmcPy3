"""
Jordyn Kim

this files pulls out the coordinates and the energy for each energy calculation for scan files and stores it in a xyz file
"""


import sys
import numpy as np
import FlexibleSPCEPotentialComplete as fpc
import matplotlib.pyplot as plt

def get_SCF(outputFile):
        readFile = open(outputFile,'r')

        SCFlist = []

        for line in readFile:
                if 'SCF Done' in line:
                        line = line.split()
                        SCFlist.append(float(line[4]))

        readFile.close()
        return np.array(SCFlist)


def get_coords2(outputFile, SCFlist, numbatoms):
        readFile = open(outputFile,'r')
        elementLookUp = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'X']
        countdown = 0
        scfcount = 0 


        # for water 1
        coordsWater1O = []
        coordsWater1H = []
        coordsWater1H2 = []
        # for water 2
        coordsWater2O = []
        coordsWater2H = []
        coordsWater2H2 = []



        for line in readFile:
                if "Standard orientation:" in line:
                        countdown = 6 + numbatoms
                        # print("New Run")
                        # print(SCFlist[scfcount])
                        scfcount += 1
                        atomcount = 0
                if countdown > numbatoms:
                        countdown = countdown - 1
                if countdown <= numbatoms and countdown > 0:
                        linelist = line.split()
                        x = float(linelist[3])*1.88973 
                        y = float(linelist[4])*1.88973 
                        z = float(linelist[5])*1.88973 

                        if countdown == 6:
                                coordsWater1O.append([x, y, z])
                        elif countdown == 5:
                                coordsWater1H.append([x, y, z])
                        elif countdown == 4:
                                coordsWater1H2.append([x, y, z])
                        elif countdown == 3:
                                coordsWater2O.append([x, y, z])
                        elif countdown == 2:
                                coordsWater2H.append([x, y, z])
                        else:
                                coordsWater2H2.append([x, y, z])
                        countdown = countdown - 1

        readFile.close()
        return coordsWater1O, coordsWater1H, coordsWater1H2, coordsWater2O, coordsWater2H, coordsWater2H2

def main(fileName, numbatoms):

        SCFlist = get_SCF(fileName)
        SCFlist = SCFlist - np.min(SCFlist)
        # get_coords(fileName, fileName2, SCFlist)
        coordsWater1O, coordsWater1H, coordsWater1H2, coordsWater2O, coordsWater2H, coordsWater2H2 = get_coords2(fileName, SCFlist, numbatoms)

        VinterSumList = []
        coloumbicEnergySumList = []
        lennardJonesSumList = []
        xaxis = []

        for i in range(len(coordsWater1O)):
                # print(coordsWater1O[i])
                atom1 = [coordsWater1O[i], coordsWater1H[i], coordsWater1H2[i]]

                atom2 = [coordsWater2O[i], coordsWater2H[i], coordsWater2H2[i]]

                atom1 = np.array(atom1)
                atom2 = np.array(atom2)

                VinterSum, coloumbicEnergySum, lennardJonesSum = fpc.PotentialEnergyTwoWaters(atom1, atom2)
                VinterSumList.append(VinterSum)
                coloumbicEnergySumList.append(coloumbicEnergySum)
                lennardJonesSumList.append(lennardJonesSum)
                xaxis.append(i)

        return xaxis, VinterSumList, coloumbicEnergySumList, lennardJonesSumList

if __name__ == '__main__':
        fileName = sys.argv[1]
        # fileName2 = sys.argv[2]
        numbatoms = input("How many atoms are in your input besides the dummy atom? ")
        numbatoms = int(numbatoms)
        xaxis, VinterSumList, coloumbicEnergySumList, lennardJonesSumList = main(fileName, numbatoms)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.scatter(xaxis, VinterSumList)
        ax1.set_xlabel(F'Run (Angstrom)')
        ax1.set_ylabel('Energy (a.u.)')
        ax1.set_title(F'{fileName}')
        plt.show()
