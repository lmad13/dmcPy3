"""
try to find functional fit of .log files with 2D scan variables
fixed the wireframe by using plot_trisurf
combineed the reflection of the graph over the other side of the angle
added a3*np.sin(b3*(y-h3)) + a4*np.sin(b4*(y-h4)) + a3*np.sin(b4*(y-h4))
repeat the pattern on the ends, to make a better fit for the data

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
    
    editedyreverseleft = -1*np.flip(editedy)
    
    #middle two reflections that truly represent the PES
    combinex = np.concatenate((editedx, editedxreverse[19:]))  
    combiney = np.concatenate((editedy, editedyreverse[19:]))
    combinez = np.concatenate((editedz, editedzreverse[19:]))
    
    #combined with left reflection (below 0 rads) that does not truly represent the PES
    leftflipx = np.flip(combinex)
    leftflipy = -1* np.flip(combiney) 
    leftflipz = np.flip(combinez)
    
    combinexleft = np.concatenate((leftflipx[:-19], combinex))
    combineyleft = np.concatenate((leftflipy[:-19], combiney))
    combinezleft = np.concatenate((leftflipz[:-19], combinez))
    
    #combined with right reflection (above 2pi rads) that does not truly represent the PES
    rightflipx = np.flip(combinex)
    rightflipy = 4*np.pi - np.flip(combiney) 
    rightflipz = np.flip(combinez)    

    combinexright = np.concatenate((combinexleft, rightflipx[19:]))
    combineyright = np.concatenate((combineyleft, rightflipy[19:]))
    combinezright = np.concatenate((combinezleft, rightflipz[19:]))
    

    def model(var, w, a0, a1, a2, a3, a4, a5, a6,a7, a8, a9, c1, h):
        x, y = var
        gy = a0 + a1*np.cos(w*y) + a2*np.cos(2*w*y) + a3*np.cos(3*w*y) + a4*np.cos(4*w*y) + a5*np.cos(5*w*y) +  a6*np.cos(6*w*y) +  a7*np.cos(7*w*y) +  a8*np.cos(8*w*y) + a9*np.cos(9*w*y)
        fx = np.exp(-1*c1*x)
        return gy* fx + h
    
    initial_fit = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0 , 1.0]
    
    fit = curve_fit(model, (combinexright, combineyright), combinezright, initial_fit)
        
    ans, cov = fit
    
    fit_w, fit_a0, fit_a1, fit_a2, fit_a3, fit_a4, fit_a5, fit_a6, fit_a7, fit_a8, fit_a9, fit_c1, fit_h = ans
    
    #This is on the entire reflected PES that don't necessarily need it
    #modeledz =  model((combinexright, combineyright), fit_w, fit_a0, fit_a1, fit_a2, fit_a3, fit_a4, fit_a5, fit_a6, fit_a7, fit_a8, fit_c1, fit_h) 
    
    """
    residuals = combinezright - modeledz
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((combinezright - np.mean(combinezright))**2)
    r_squared = 1- (ss_res/ss_tot)
    """
    
    
    #This is on just the PES that need it
    modeledz =  model((combinex, combiney), fit_w, fit_a0, fit_a1, fit_a2, fit_a3, fit_a4, fit_a5, fit_a6, fit_a7, fit_a8, fit_a9, fit_c1, fit_h)
    

    residuals = combinez - modeledz
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((combinez - np.mean(combinez))**2)
    r_squared = 1- (ss_res/ss_tot)

    
    
        
    print('fit_w: ', fit_w)
    print('fit_a0: ',fit_a0)
    print('fit_a1: ', fit_a1)
    print('fit_a2: ', fit_a2)
    print('fit_a3: ', fit_a3)
    print('fit_a4: ', fit_a4)
    print('fit_a5: ', fit_a5)
    print('fit_a6: ', fit_a6)
    print('fit_a7: ', fit_a7)
    print('fit_a8: ', fit_a8)
    print('fit_a9: ', fit_a9)
    print('fit_c1: ', fit_c1)
    print('fit_h: ', fit_h)
    print('SS_res:', ss_res)
    print('SS_tot:', ss_tot)
    print('r_squared:', r_squared)
    
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    
    #This is on the entire reflected PES that don't necessarily need it
    #ax1.scatter(combinexright, combineyright, combinezright, c=combinezright, label="Original Data", cmap='RdBu')
    #ax1.plot_trisurf(combinexright, combineyright, modeledz, color='black', alpha=0.3)
    
    #This is on just the PES that need it
    ax1.scatter(combinex, combiney, combinez, c=combinez, label="Original Data", cmap='RdBu')
    ax1.plot_trisurf(combinex, combiney, modeledz, color='black', alpha=0.3)    
    
    ax1.set_xlabel(F'{var1in} (Angstrom)')
    ax1.set_ylabel(F'{var2in} (Radians)')
    ax1.set_zlabel('Energy (a.u.)')
    ax1.set_title(F'{fileName}')
    plt.show()
