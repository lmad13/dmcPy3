import dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys


#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)
au2wn=219474.63


#The number of walkers

nWalkersmin=1900
nWalkersmax=2000

walkers=np.array([i+1 for i in range(nWalkersmin, nWalkersmax)])
vavg=np.zeros(nWalkersmax-nWalkersmin)

count=0

for i in range(nWalkersmin,nWalkersmax):
    H2Wfn=dmc.wavefunction(i+1,'harmonic',plotting=False)

    TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()

    print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+' cm^-1')

#Important parameters for the MC simulation
    nReps=3  #Repeat the DMC simulation 3 times
    nEquilibrationSteps=10 #Initially equilibrate the simulation with 12000 diffusion steps 
    nSteps=200  #Make 2000 Diffusion steps in each simulation

#don't worry about this just yet
    nDesSteps=75  
    nRepsDW=2000


    print('important parameters:')
    print('Equilibrating for '+str(nEquilibrationSteps))
    print('The number of diffusion steps in the simulation is '+str(nSteps))
    print('The time step is ' +str(H2Wfn.dtau))


    H2Wfn.setX(H2Wfn.xcoords+.2)


#first equilibrate:
    v,pop,x0Equilibration,d=H2Wfn.propagate(H2Wfn.xcoords,nEquilibrationSteps)

    vcombined = np.zeros(nReps)


#then simulate n times:
    print('nWalkers: ', i+1)
    for n in range(nReps):
    #propagate

        
        vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps,plotWalkers=False,printFlag=False)
        vref_0=np.array(vref_0)
    

        #print("Repetition: " + str(n))
        #print("Average Energy:")
        #print(np.average(vref_0))
        vcombined[n]=np.average(vref_0)
    

    vavg[count]= np.average(vcombined)
    count=count+1

print('vavg: ',vavg)
print('walkers: ', walkers)

plt.figure()
plt.plot(walkers, vavg, 'o')
plt.xlabel('Number of Walkers')
plt.ylabel('V average')
plt.title('V average vs # of Walkers')
plt.show()

