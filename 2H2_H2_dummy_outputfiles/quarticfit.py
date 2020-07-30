#Quartic Potential Energy Function Fit

"""
graph SCF energies for 2H2_H2_dummy_input_0.log and 2H2_H2_dummy_input_1.log
"""

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import sys
from scipy.optimize import curve_fit

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
            if len(listline) == 3:
                if len(listline[0]) < 6 & len(listline[1]) > 2:
                    b6 = float(listline[1])
                    variablelist.append(b6)
    
    readFile.close()
    return np.array(variablelist)

x = PrintArrayb6(fileName)
y = PrintArraySCF(fileName)

editedx = x - x[0]
editedy = y - np.min(y)

editedxright = -1*editedx

editedxleft = np.flip(editedx)
editedyleft = np.flip(editedy)

combinex = np.concatenate((editedxleft, editedxright[1:]))
combiney = np.concatenate((editedyleft, editedy[1:]))

print("B6: ", combinex)
print("SCF energy: ", combiney)


def qmodel(x, qa, qb, qc, qd, qe):
    return qa*(x**4)+qb*(x**3)+qc*(x**2)+qd*(x)+qe

qfit = curve_fit(qmodel, combinex, combiney)

qans, qcov = qfit
fit_qa, fit_qb, fit_qc, fit_qd, fit_qe = qans


modeledy = qmodel(combinex, fit_qa, fit_qb, fit_qc, fit_qd, fit_qe)

residuals = combiney - modeledy
ss_res = np.sum(residuals**2)
ss_tot = np.sum((combiney - np.mean(combiney))**2)
r_squared = 1- (ss_res/ss_tot)

print('---Fit Parameters---')

print('qa: ', fit_qa)
print('qb: ', fit_qb)
print('qc: ', fit_qc)
print('qd: ', fit_qd)
print('qe: ', fit_qe)

print('SS_res:', ss_res)
print('SS_tot:', ss_tot)
print('r_squared:', r_squared)


au2wn = 219474.64

plt.figure()
plt.plot(combinex*0.529177249, combiney,'o',label='Original Data')
plt.plot(combinex*0.529177249, modeledy, '--', label='Quartic Fit')
plt.xlabel('Dummy-Oxygen Bond Length (Angstrom)')
plt.ylabel('SCF Energy (Hartree)')
plt.legend()
plt.title('Energy vs Bond Length')
plt.show()
