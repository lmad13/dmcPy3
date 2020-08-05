"""
This code is a copy of distancematrix.py so I can use each element in the matrix as paramters for a new function

This code only works for 2 H2O and 1 H2 since I use brute force to pull out elements

"""

import numpy as np
import matplotlib.pyplot as plt
import sys

print("The arguments are: " + str(sys.argv))
if len(sys.argv) < 2:
    print('Uh Oh! You did not provide a command line argument that had the name of the file!')
    exit()
else:
    fileName = sys.argv[1]
    print("reading " + fileName)
    
numbatoms = input("How many atoms are in your input besides the dummy atom? ")
numbatoms = int(numbatoms)
    
def get_coords(outputFile):
    readFile = open(outputFile,'r')
    elementLookUp = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'X']
    countdown = 0
    scfcount = 0  
    runcount = 0
    coordscol = []


    for line in readFile:
        if "Standard orientation:" in line:
            countdown = 6 + numbatoms
            scfcount += 1
            atomcount = 0
            runcount +=1
            coordsrow = [] 
        if countdown > numbatoms:
            countdown = countdown - 1
        if countdown <= numbatoms and countdown > 0:
            linelist = line.split()
            xyzcoords = [float(linelist[3]),float(linelist[4]),float(linelist[5])]
            coordsrow.append(xyzcoords)
            countdown = countdown - 1
            if countdown == 1:
                coordscol.append(coordsrow)
                        
                        

    readFile.close()
    return np.array(coordscol)

coords = get_coords(fileName)
num_rows, num_cols, num_ind = coords.shape

def get_rijmatrix(coordmatrix, numb_rows, numb_atoms):
    coordmatrix = coordmatrix
    numbruns = numb_rows
    numbatoms = numb_atoms
    rijmatrixset = []
    for run in range(numbruns):
        print("Scan number ", run)
        emptyrijmatrix = [[0.0]*numbatoms for i in range(numbatoms)]
        emptyrijmatrix = np.array(emptyrijmatrix)
        for i in range(numbatoms):
            for j in range(numbatoms):
                ix = coordmatrix[run, i, 0]
                iy = coordmatrix[run, i, 1]
                iz = coordmatrix[run, i, 2]
                jx = coordmatrix[run, j, 0]
                jy = coordmatrix[run, j, 1]
                jz = coordmatrix[run, j, 2]

                distance = np.sqrt((ix-jx)**2+(iy-jy)**2+(iz-jz)**2)
                emptyrijmatrix[i,j]=distance
        print(emptyrijmatrix)
        rijmatrixset.append(emptyrijmatrix)
        np.savetxt("distancematrix"+str(run)+".txt",emptyrijmatrix)
                
    return np.array(rijmatrixset)
          
rijmatrix = get_rijmatrix(coords, num_rows, numbatoms)
rij_runs, rij_cols, rij_rows = rijmatrix.shape

def select_rijelements(matrix, numbruns, numbatoms):
    matrix = matrix
    numbruns = numbruns
    numbatoms = numbatoms
    O1 = np.array([0])
    HH1 = np.array([1,2])
    O2 = np.array([3])
    HH2 = np.array([4,5])
    HH3 = np.array([6,7])
    for run in range(numbruns):
        print("Scan number ", run)
        vfunct = []
        for i in range(numbatoms):
            for j in range(numbatoms):
                if i < j:
                    if i in O1:
                        if j in O2:
                            vfunct.append('a')
                        elif j in HH2:
                            vfunct.append('b')
                        elif j in HH3:
                            vfunct.append('c')
                    elif i in HH1:
                        if j in O2:
                            vfunct.append('d')
                        elif j in HH2:
                            vfunct.append('e')
                        elif j in HH3:
                            vfunct.append('f')
                    elif i in O2:
                        if j in HH3:
                            vfunct.append('g')
                    elif i in HH2:
                        if j in HH3:
                            vfunct.append('h')
        print(np.array(vfunct))
    return
                            

select_rijelements(rijmatrix, rij_runs, rij_cols)
    



            