import  dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys



#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)
au2wn=219474.63


#The number of walkers
nWalkers=2000

H2H2OWfn=dmc.wavefunction(nWalkers,'morse',molecule='H2H2O')

#TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()

#print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+' cm^-1')

#Important parameters for the MC simulation
nReps=3  #Repeat the DMC simulation 3 times
nEquilibrationSteps=1000 #Initially equilibrate the simulation with 12000 diffusion steps 
nSteps=2000  #Make 2000 Diffusion steps in each simulation

nDesSteps=75  
nRepsDW=5

print('important parameters:')

print('Equilibrating for '+str(nEquilibrationSteps))
print('The number of diffusion steps in the simulation is '+str(nSteps))
print('The time step is ' +str(H2H2OWfn.dtau))


H2H2OWfn.setX(H2H2OWfn.xcoords+.2)

#first equilibrate:
v,pop,x0Equilibration,d=H2H2OWfn.propagate(H2H2OWfn.xcoords,nEquilibrationSteps)

#then simulate n times:
for n in range(nReps):
    #propagate
    vref_0,pop,x0,d=H2H2OWfn.propagate(x0Equilibration,nSteps)
    vref_0=np.array(vref_0)
    print("Repetition: " + str(n))
    print("Average Energy:")
    print(np.average(vref_0))
    print(str(np.average(vref_0)*(219474.63))+" cm-1")
    

    
