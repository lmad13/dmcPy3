import dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys


#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)
au2wn=219474.63


#The number of walkers

nWalkersmin=500
nWalkersmax=2000

walkers=np.array([i+1 for i in range(nWalkersmin, nWalkersmax)])

correctwalkers=walkers[walkers%50==0]
vavg=np.zeros(correctwalkers.size)
stdvavg=np.zeros(correctwalkers.size)

count=0

for i in range(correctwalkers.size):
    thiswalker = correctwalkers[i]
    H2Wfn=dmc.wavefunction(thiswalker,'harmonic',plotting=False)

    TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()

    print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+' cm^-1')

#Important parameters for the MC simulation
    nReps=10  #Repeat the DMC simulation 3 times
    nEquilibrationSteps=10 #Initially equilibrate the simulation with 12000 diffusion steps 
    nSteps=200  #Make 2000 Diffusion steps in each simulation

#don't worry about this just yet
    nDesSteps=75  
    nRepsDW=2000


    print('important parameters:')
    print('Equilibrating for '+str(nEquilibrationSteps))
    print('The number of diffusion steps in the simulation is '+str(nSteps))
    print('The time step is ' +str(H2Wfn.dtau))
    print('number of reps is ', nReps)

    H2Wfn.setX(H2Wfn.xcoords+.2)


#first equilibrate:
    v,pop,x0Equilibration,d=H2Wfn.propagate(H2Wfn.xcoords,nEquilibrationSteps)

    vcombined = np.zeros(nReps)


#then simulate n times:
    print('nWalkers: ', thiswalker)
    for n in range(nReps):
    #propagate

        
        vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps,plotWalkers=False,printFlag=False)
        vref_0=np.array(vref_0)
    

        #print("Repetition: " + str(n))
        #print("Average Energy:")
        #print(np.average(vref_0))
        vcombined[n]=np.average(vref_0)
        

    vavg[count]= np.average(vcombined)
    stdvavg[count]=np.std(vcombined)
    count=count+1

print('vavg: ',vavg)
print('walkers: ', correctwalkers)
print('std', stdvavg)

plt.figure()
plt.errorbar(correctwalkers, vavg, yerr=stdvavg, color='orange', ecolor='green', capsize = 5, capthick=2, fmt='o')
#plt.plot(correctwalkers, vavg, 'o')
plt.xlabel('Number of Walkers')
plt.ylabel('V average')
plt.title('V average vs # of Walkers')

#plt.figure()
#plt.plot(correctwalkers, stdvavg, 'o')
#plt.xlabel('Number of Walkers')
#plt.ylabel('Std V average')
#plt.title('Std V average vs # of Walkers')


plt.show()

