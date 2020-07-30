import sys
import numpy as np
import matplotlib.pyplot as plt

#graphSCFvsB6.py 2H2_H2_dummy_input.log
#graph SCF vs B6
#obtain 2 log files to combine left(approaching) and right(moving away) scan                                                                              #if one file, then doesn't work   
print("The arguments are: " + str(sys.argv))
if len(sys.argv) < 2:
    print('Uh Oh! You did not provide a command line argument that had the name of the file!')
    exit()
else:
    fileName = sys.argv[1]
    print("reading " + fileName)

#obtain index of SCF by collecting a list of index before optimized parameters 
def PrintArraySCF(fileName):
    readFile = open(fileName, 'r')
    printThisLine = False
    scflist = []
    count = -1 #start with -1 because I am getting index that has optimized paramters but want SCF, so subtract 1 from the index to get where SCF energy index is
    optcount = [] #initialize an empty array to collect where SCF energy is

    for line in readFile:
        if 'SCF Done' in line:
            printThisLine = True
        if 'NFock=' in line:
            printThisLine = False
        if 'Optimized Parameters' in line:
            printThisLine = True
        if 'Angstroms and Degrees' in line:
            printThisLine = False
        if printThisLine:
            if 'Optimized Parameters' in line:
                optcount.append(count)
            count+=1


    return optcount

#obtain array of where scf is
scfcount = PrintArraySCF(fileName)

#This obtains all the values of SCF's 
def PrintRealSCF(fileName, scfcount):#scfcount obtains where scf is  
    readFile = open(fileName, 'r')
    printThisLine = False
    scflist = []
    count2 = 0 #will use this to run through the index 
    optcount = scfcount

    for line in readFile:
        if 'SCF Done' in line:
            printThisLine = True
        if 'NFock=' in line:
            printThisLine = False
        if 'Optimized Parameters' in line:
            printThisLine = True
        if 'Angstroms and Degrees' in line:
            printThisLine = False
        if printThisLine:
            if count2 in optcount:#if the index(with SCF's and optimized) is in the array that only gives out scfs  
                correct = line.split()
                scflist.append(float(correct[4]))
            count2=1+count2#next index

    readFile.close()
    return np.array(scflist)


def PrintBondLength(fileName):
    readFile = open(fileName, 'r')
    printThisLine = False
    blist = []
    for line in readFile:
        if ' !       B6' in line:
            if 'Scan' not in line:
                printThisLine = True
        if ' !       A6' in line:
            printThisLine = False
        if printThisLine:
            lengthline = line.split()
            blist.append(float(lengthline[2]))

    readFile.close()
    return blist


#length in angstrom
b6 = np.array(PrintBondLength(fileName))
#length in bohr
editedb6=b6*1.8897259886
#original energy in Hartree
scf = PrintRealSCF(fileName, scfcount)
lowestscf = np.amin(scf)
#set lowest as 0 energy
editedscf = scf-lowestscf

print('---SCF Done---')
print(scf[0:-1])
print(editedscf[0:-1])
print('Lowest energy is: ', lowestscf)


print('---B6---')
print(editedb6[0:-1])



plt.figure()
plt.plot(editedb6[0:-1], editedscf[0:-1],'o')
plt.xlabel('B6 (Bohr)')
plt.ylabel('SCF Energy (Hartree)')
plt.title('Energy vs Bond Length')
plt.show()
