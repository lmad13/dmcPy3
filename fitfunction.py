"""
try to find functional fit of .log files with 2D scan variables
fixed the wireframe by using plot_trisurf
combineed the reflection of the graph over the other side of the angle
added a3*np.sin(b3*(y-h3))
works for 2H2_H2_dummy_input_5.log
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import sys
from scipy import interpolate
import pylab as py
from scipy.optimize import curve_fit
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


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

    editedx = x - x[0]
    editedy = np.radians(y[0]) - np.radians(y)  #change angles in angles to radians
    editedz = z - np.amin(z)
    
    editedxreverse = np.flip(editedx)
    editedyreverse = 2*np.pi - np.flip(editedy)
    editedzreverse = np.flip(editedz)
    

    combinex = np.concatenate((editedx, editedxreverse[20:]))
    combiney = np.concatenate((editedy, editedyreverse[20:]))
    combinez = np.concatenate((editedz, editedzreverse[20:]))

    
    def model(variables, c1, a2, b2, h2, a3, b3, h3, k):
        x, y = variables
        gy = a2*np.sin(b2*(y-h2)) + a3*np.sin(b3*(y-h3))
        fx = np.exp(-1*c1*x)
        return fx * gy + k
    
    initial_fit = [3.2682,0.0047,3.0619,-0.4494,0.0203,0.4708,-13.5383,-0.0044]
    
    fit = curve_fit(model, (combinex, combiney), combinez, initial_fit)
    
    ans, cov = fit
    
    fit_c1, fit_a2, fit_b2, fit_h2, fit_a3, fit_b3, fit_h3, fit_k = ans
    
    modeledz =  model((combinex, combiney), fit_c1, fit_a2, fit_b2, fit_h2, fit_a3, fit_b3, fit_h3, fit_k)
    
    residuals = combinez - modeledz
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((combinez - np.mean(combinez))**2)
    r_squared = 1- (ss_res/ss_tot)
    
    print('c1: ', fit_c1)
    print('a2: ', fit_a2)
    print('b2: ', fit_b2)
    print('h2: ', fit_h2)
    print('a3: ', fit_a3)
    print('b3: ', fit_b3)
    print('h3: ', fit_h3)
    print('k: ', fit_k)
    print('SS_res:', ss_res)
    print('SS_tot:', ss_tot)
    print('r_squared:', r_squared)
    
 
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    ax1.scatter(combinex, combiney, combinez, c=combinez, cmap='RdBu') #original data
    ax1.plot_trisurf(combinex, combiney, modeledz, color='black', alpha=0.2)
    ax1.set_xlabel(F'{var1in} (Angstrom)')
    ax1.set_ylabel(F'{var2in} (Radians)')
    ax1.set_zlabel('Energy (a.u.)')
    ax1.set_title(F'{fileName}')
    plt.show()

        