import dmc1D2 as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys

#graph v average vs. dtau with standard deviations

#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)                              
au2wn=219474.63


#The number of walkers                                                                                        
nWalkers=2000
rangedtau = 5

#V averages rangedtau many times                                                                              
vaverages = np.zeros(rangedtau)

#std V rangedtau many times                                                                                   
stdvaverages = np.zeros(rangedtau)

#dtaus                                                                                                        
dtaus = np.array([i+1 for i in range(rangedtau)])


for i in range(rangedtau):
    H2Wfn=dmc.wavefunction(nWalkers,'harmonic',plotting=False, ndtau=i+1)

    TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()
    print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+' cm^-1')

#Important parameters for the MC simulation                                                                   
    nReps=3  #Repeat the DMC simulation 3 times                                                               
    nEquilibrationSteps=10 #Initially equilibrate the simulation with 12000 diffusion steps                   
    nSteps=int(200/(i+1))  #Make 2000 Diffusion steps in each simulation                                      
    nReps_vaverages = np.zeros(nReps) #make V averages nReps many times                                       
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
    v,pop,x0Equilibration,d=H2Wfn.propagate(H2Wfn.xcoords,nEquilibrationSteps)

#then simulate n times:                                                                                       


    for n in range(nReps):

    #propagate                                                                                                
        vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps,plotWalkers=False,printFlag=False)
        vref_0=np.array(vref_0)

        average_energy = np.average(vref_0)
        standard_deviation = np.std(vref_0)


        print("Repetition: " + str(n))
        print("Average Energy:")
        print(average_energy)
        print('Standard Deviation of Average Energy:')
        print(standard_deviation)
        print(' ')

        nReps_vaverages[n] = average_energy

    print('Average of Average Energy for dtau {}:'.format(i+1))
    print(np.average(nReps_vaverages))
    vaverages[i]=np.average(nReps_vaverages)
    stdvaverages[i]=(np.std(nReps_vaverages))
    print(' \n-------------------------\n ')


print('List of dtaus:')
print(dtaus)
print('List of Average of Average Energy')
print(vaverages)
print('List of Std of Average Energy')
print(stdvaverages)


plt.errorbar(dtaus, vaverages, yerr=stdvaverages, color='orange', ecolor='green', capsize = 5, capthick=2, fmt='o')
plt.xlabel('time step, dtau')
plt.ylabel('V average')
plt.title('V average vs dtau')
plt.xticks(dtaus)
plt.show()

