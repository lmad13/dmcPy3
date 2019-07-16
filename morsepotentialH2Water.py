import dmc1Dmorse as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys



#Conversion factor of atomic units of energy to wavenumber (inverse centimeters)
au2wn=219474.63


#The number of walkers
nWalkers=10000

H2Wfn=dmc.wavefunction(nWalkers,'Morse',plotting=False)

TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()

print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+' cm^-1')

#Important parameters for the MC simulation
nReps=3  #Repeat the DMC simulation 3 times
nEquilibrationSteps=1200 #Initially equilibrate the simulation with 12000 diffusion steps 
nSteps=2000  #Make 2000 Diffusion steps in each simulation

#don't worry about this just yet
nDesSteps=75  
nRepsDW=2000

print('important parameters:')

print('Equilibrating for '+str(nEquilibrationSteps))
print('The number of diffusion steps in the simulation is '+str(nSteps))
print('The time step is ' +str(H2Wfn.dtau))


H2Wfn.setX(H2Wfn.xcoords+6.2)

#first equilibrate:
v,pop,x0Equilibration,d=H2Wfn.propagate(H2Wfn.xcoords,nEquilibrationSteps,plotWalkers=False,printFlag=False)

vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps)
vref_0=np.array(vref_0)
print("Repetition: 0")
print("Average Energy:")
print(np.average(vref_0))
#graph normalized histogram of x0

hist, bin_edges = np.histogram(x0, bins = 100)
norm_hist = hist/x0.size
norm_bin = (bin_edges[:-1]+bin_edges[1:])/2

plt.figure(2)
plt.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run #{}'.format(1))
plt.grid()

#graph cumulative histogram of x0
cum_hist = np.cumsum(norm_hist)

plt.figure(3)
plt.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run#{}'.format(1))
plt.grid()
hiscombined = norm_hist

#graph vref vs time steps
plt.figure(4)
steps = np.array([i+1 for i in range(nSteps)])
plt.scatter(steps, vref_0, label='V ref for run {}'.format(1), s=700/nSteps)

#modify ticks                                                                                                                 
preticks = steps
mod = preticks % 100 == 0
ticks = preticks[mod]

#graph psi squared
plt.figure(1)
plt.plot(norm_bin, norm_hist**2,'o',color='{}'.format(1.0/(2**(1))), label='Run #{}'.format(1))

#then simulate n times:
for n in range(nReps-1):
    #propagate
    vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps)
    vref_0=np.array(vref_0)
    print("Repetition: " + str(n+1))
    print("Average Energy:")
    print(np.average(vref_0))
    print(x0)
#plot vref_0 over time in comparison to v avarage                                                             
    plt.figure(4)
    steps = np.array([i+1 for i in range(nSteps)])
    plt.scatter(steps, vref_0, label='V ref for run {}'.format(2+n), s=700/nSteps)
    #plot normalized histogram of population's x value                                                                                                                                                                 
    hist, bin_edges = np.histogram(x0, bins = bin_edges)
    norm_hist = hist/x0.size
    norm_bin = (bin_edges[:-1]+bin_edges[1:])/2
    plt.figure(2)
    plt.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run #{}'.format(n+2))
    plt.grid()
#plot psi squared
    plt.figure(1)
    plt.plot(norm_bin, norm_hist**2,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))
    plt.grid()

    newhiscombined = hiscombined + norm_hist
    hiscombined = newhiscombined
    
    #plot cumulative hist
    cum_hist = np.cumsum(norm_hist)
    plt.figure(3)
    plt.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))
    plt.grid()



#graph averaged histogram
avgnormedhist = hiscombined/3
plt.figure(2)
plt.plot(norm_bin,avgnormedhist, label='Avg Psi')
plt.xlabel('x (Bohr)')
plt.ylabel('frequency')
plt.grid()
plt.title('Cumulative Wavefunction for {} Runs'.format(nReps))
plt.legend(loc='best')

#graph psi squared
plt.figure(1)
plt.plot(norm_bin,avgnormedhist**2, label='Avg Psi Squared')
plt.figure(1)
plt.xlabel('x (Bohr)')
plt.ylabel('Probability')
plt.grid()
plt.title('Cumulative Wavefunction for {} Runs'.format(nReps))
plt.legend(loc='best')

#graphed cumulative histogram
plt.figure(3)
plt.xlabel('x (Bohr)')
plt.ylabel('frequency')
plt.grid()
plt.title('Wavefunction for {} Runs'.format(nReps))
plt.legend(loc='best')

#graph vref
plt.figure(4)
plt.xticks(ticks)
plt.axhline(y = np.average(vref_0), color='{}'.format(1.0/(2**(1))),label='V average for run {}: {}'.format(1,np.average(vref_0)))
plt.xlabel('Step')
plt.ylabel('V')
plt.title('V ref over time for {} Runs'.format(nReps))
plt.legend(loc='best')
plt.xlim(0,nSteps)


plt.show()
