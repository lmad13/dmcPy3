import dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys

#graph psi squared with different time steps
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

#propagate
vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps)
vref_0=np.array(vref_0)
print("Repetition: 0")
print("Average Energy:")
print(np.average(vref_0))

#graph normalized histogram of x0                                                                                                 

hist, bin_edges = np.histogram(x0, bins = 100)
norm_hist = hist/x0.size
norm_bin = (bin_edges[:-1]+bin_edges[1:])/2

#desc propagate:
descsteps = 200
initsteprange = np.array([n for n in range(descsteps)])
stepmod = initsteprange % 10 == 0
steprange = initsteprange[stepmod]

for i in steprange:
    vrefdw, popdw, xdw, descendants = H2Wfn.propagate(x0, i)
    #graph psi squared
    histsqrt, sqrtbin_edges = np.histogram(x0, bins=bin_edges, weights=descendants)
    norm_hist2 = histsqrt * norm_hist
    norm_histsqrt = norm_hist2/np.sum(norm_hist2) 
    plt.figure(1)
    plt.plot(norm_bin, norm_histsqrt,'o', label='Run #{}'.format(i))

plt.figure(1)
plt.legend(loc='best')
plt.xlabel('x (Bohr)')
plt.ylabel('Probability')
plt.grid()
plt.title('Psi Squared With Different Steps')
plt.show()
