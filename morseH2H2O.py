import dmc1Dmorse as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys

au2wn=219474.63

nWalkers=10000

H2Wfn=dmc.wavefunction(nWalkers,'Morse',plotting=False)

TheoreticalOmega0=H2Wfn.getTheoreticalOmega0()
print('the theoretical frequency for H2 vibration is: '+str(TheoreticalOmega0)+\
' cm^-1')

nReps=3
nEquilibrationSteps=1200
nSteps=2000

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

hist, bin_edges = np.histogram(x0, bins = 100)
norm_hist = hist/x0.size
norm_bin = (bin_edges[:-1]+bin_edges[1:])/2

plt.figure(1)
plt.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run #{}'.format(1))

cum_hist = np.cumsum(norm_hist)
plt.figure(2)
plt.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run#{}'.format(1))

hiscombined = norm_hist

for n in range(nReps-1):
    vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps)
    vref_0=np.array(vref_0)
    print("Repetition: " + str(n+1))
    print("Average Energy:")
    print(np.average(vref_0))
    print(x0)

    hist, bin_edges = np.histogram(x0, bins = bin_edges)
    norm_hist = hist/x0.size
    norm_bin = (bin_edges[:-1]+bin_edges[1:])/2
    plt.figure(1)
    plt.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run #{}'.format(n+2))

    newhiscombined = hiscombined + norm_hist
    hiscombined = newhiscombined

    cum_hist = np.cumsum(norm_hist)
    plt.figure(2)
    plt.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))

avgnormedhist = hiscombined/3
plt.figure(1)
plt.plot(norm_bin,avgnormedhist)

plt.figure(1)
plt.xlabel('x (Bohr)')
plt.ylabel('frequency')
plt.grid()
plt.title('Cumulative Wavefunction for {} Runs'.format(nReps))
plt.legend(loc='best')

plt.figure(2)
plt.xlabel('x (Bohr)')
plt.ylabel('frequency')
plt.grid()
plt.title('Wavefunction for {} Runs'.format(nReps))
plt.legend(loc='best')


plt.show()

