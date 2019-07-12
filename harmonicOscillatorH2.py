import dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys


#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)
au2wn=219474.63


#The number of walkers
nWalkers=2000


H2Wfn=dmc.wavefunction(nWalkers,'harmonic',plotting=False)

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

x0combined = []
bincombined = np.empty(0)


#then simulate n times:
for n in range(nReps):
    #propagate                                                                                        
    vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps,plotWalkers=False,printFlag=False)
    vref_0=np.array(vref_0)
    

    print("Repetition: " + str(n))
    print("Average Energy:")
    print(np.average(vref_0))
#plot vref_0 over time in comparison to v avarage                              
                                                                                    
    plt.figure(1)
    steps = np.array([i+1 for i in range(nSteps)])
    steps = np.array([i+1 for i in range(nSteps)])
    plt.scatter(steps, vref_0, label='V ref for run {}'.format(n+1), s=700/nSteps)

    #modify ticks                                                                   
    preticks = steps 
    mod = preticks % 100 == 0
    ticks = preticks[mod]

    plt.xticks(ticks)
    plt.axhline(y = np.average(vref_0), color='{}'.format(1.0/(2**(n+1))),label='V average for run {}: {}'.format(n+1,np.average(vref_0)))
    plt.xlabel('Step')
    plt.ylabel('V')
    plt.title('V ref over time for {} Runs'.format(nReps))
    plt.legend(loc='best')
    plt.xlim(0,nSteps)
    
    #plot normalized histogram of population's x value                                                          
    hist, bin_edges = np.histogram(x0, bins = 100, range=(-1.5,1.5))
    norm_hist = hist/x0.size
    

    plt.figure(2)
    plt.plot((bin_edges[:-1]+bin_edges[1:])/2, norm_hist,'o',color='{}'.format(1.0/(2**(n+1))), label='Run #{}'.format(n+1))
    plt.xlabel('x (Bohr)')
    plt.ylabel('frequency')
    plt.grid()
    plt.title('Wavefunction for {} Runs'.format(nReps))
    plt.legend(loc='best')

    
    #cumulative histogram
    #cum_hist = np.cumsum(norm_hist)
    cum_hist = np.cumsum(norm_hist)
    plt.figure(3)
    plt.plot(bin_edges[:-1], cum_hist,'o',color='{}'.format(1.0/(2**(n+1))), label='Run#{}'.format(n+1))
    plt.xlabel('x (Bohr)')
    plt.ylabel('frequency')
    plt.grid()
    plt.title('Cumulative Wavefunction for {} Runs'.format(nReps))
    plt.legend(loc='best')


    #combine x0 to whole list
    x0combined.append(norm_hist)
    #bincombined = np.concatenate((bincombined, bin_edges))
x0combined=np.array(x0combined)
avgnormedhist = np.average(x0combined,axis=0)
plt.figure(2)
plt.plot((bin_edges[:-1]+bin_edges[1:])/2,avgnormedhist)

plt.show()

