import dmc1D as dmc
import numpy as np
import matplotlib.pyplot as plt
import sys

#test of Morse

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

#initialize figures
fig1 = plt.figure(1)
fig2 = plt.figure(2)

#desc propagate:
vrefdw, popdw, xdw, descendants = H2Wfn.propagate(x0, 10)


#graph normalized histogram of x0

hist, bin_edges = np.histogram(x0, bins = 100)
norm_hist = hist/x0.size
norm_bin = (bin_edges[:-1]+bin_edges[1:])/2

ax1 = fig1.add_subplot(221)
ax1.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run #{}'.format(1))
hiscombined = norm_hist

#graph integrated histogram of x0
cum_hist = np.cumsum(norm_hist)

bx1 = fig2.add_subplot(121)
bx1.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(1))), label='Run#{}'.format(1))

#graph vref vs time steps
plt.figure(3)
steps = np.array([i+1 for i in range(nSteps)])
plt.scatter(steps, vref_0, label='V ref for run {}'.format(1), s=700/nSteps)

#modify ticks                                                                                                                 
preticks = steps
mod = preticks % 100 == 0
ticks = preticks[mod]

#graph psi squared
histsqrt, sqrtbin_edges = np.histogram(x0, bins=bin_edges, weights=descendants)
norm_hist2 = histsqrt * norm_hist
norm_histsqrt = norm_hist2/np.sum(norm_hist2) 

ax2 = fig1.add_subplot(222)
ax2.plot(norm_bin, norm_histsqrt,'o',color='{}'.format(1.0/(2**(1))), label='Run #{}'.format(1))
histsqrtcombined = norm_histsqrt

#plot cumulative histsqrt                                                                                                     
cum_histsqrt = np.cumsum(norm_histsqrt)
bx2 = fig2.add_subplot(122)
bx2.plot(norm_bin, cum_histsqrt,'o',color='{}'.format(1.0/(2**(1))), label='Run#{}'.format(1))

#then simulate n times:
for n in range(nReps-1):
    #propagate
    vref_0,pop,x0,d=H2Wfn.propagate(x0Equilibration,nSteps)
    vref_0=np.array(vref_0)
    print("Repetition: " + str(n+1))
    print("Average Energy:")
    print(np.average(vref_0))

    #descendant propagate:
    vrefdw, popdw, xdw, descendants = H2Wfn.propagate(x0, 10)

#plot vref_0 over time in comparison to v avarage                                                             
    plt.figure(3)
    steps = np.array([i+1 for i in range(nSteps)])
    plt.scatter(steps, vref_0, label='V ref for run {}'.format(2+n), s=700/nSteps)
    #plot normalized histogram of population's x value                                                                                                                                                                 
    hist, bin_edges = np.histogram(x0, bins = bin_edges)
    norm_hist = hist/x0.size
    norm_bin = (bin_edges[:-1]+bin_edges[1:])/2
    
    ax1.plot(norm_bin, norm_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run #{}'.format(n+2))

    newhiscombined = hiscombined + norm_hist
    hiscombined = newhiscombined

#plot psi squared
    histsqrt, sqrtbin_edges = np.histogram(x0, bins=bin_edges, weights=descendants)
    norm_hist2 = histsqrt *norm_hist
    norm_histsqrt = norm_hist2/np.sum(norm_hist2)

    ax2.plot(norm_bin, norm_histsqrt,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))

    newhistsqrtcombined = histsqrtcombined + norm_histsqrt
    histsqrtcombined = newhistsqrtcombined
    
    #plot cumulative hist
    cum_hist = np.cumsum(norm_hist)
    bx1.plot(norm_bin, cum_hist,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))

    #plot cumulative histsqrt
    cum_histsqrt = np.cumsum(norm_histsqrt)
    bx2.plot(norm_bin, cum_histsqrt,'o',color='{}'.format(1.0/(2**(n+2))), label='Run#{}'.format(n+2))

#graph averaged histogram
avgnormedhist = hiscombined/nReps
ax1.plot(norm_bin,avgnormedhist, label='Avg Psi')
ax1.set_xlabel('x (Bohr)')
ax1.set_ylabel('frequency')
ax1.grid()
ax1.set_title('Cumulative Wavefunction for {} Runs'.format(nReps))
ax1.legend(loc='best')

#graph psi squared
avgnormedhistsqrt = histsqrtcombined/nReps
ax2.plot(norm_bin,avgnormedhistsqrt, label='Avg Psi Squared')
ax2.set_xlabel('x (Bohr)')
ax2.set_ylabel('Probability')
ax2.grid()
ax2.set_title('Cumulative Probability for {} Runs'.format(nReps))
ax2.legend(loc='best')

#graphed cumulative histogram
bx1.set_xlabel('x (Bohr)')
bx1.set_ylabel('frequency')
bx1.grid()
bx1.set_title('Integrate Psi for {} Runs'.format(nReps))
bx1.legend(loc='best')

#graph cumulative histsqrt
bx2.set_xlabel('x (Bohr)')
bx2.set_ylabel('Probability')
bx2.grid()
bx2.set_title('Integrate Psi Squared for {} Runs'.format(nReps))
bx2.legend(loc='best')


#graph vref
plt.figure(3)
plt.xticks(ticks)
plt.axhline(y = np.average(vref_0), color='{}'.format(1.0/(2**(1))),label='V average for run {}: {}'.format(1,np.average(vref_0)))
plt.xlabel('Step')
plt.ylabel('V')
plt.title('V ref over time for {} Runs'.format(nReps))
plt.legend(loc='best')
plt.xlim(0,nSteps)

#graph two avg's into one plot
ax3 = fig1.add_subplot(212)
ax3.plot(norm_bin,avgnormedhist, label='Avg Psi')
ax3.plot(norm_bin,avgnormedhistsqrt, label='Avg Psi Squared')
ax3.set_xlabel('x (Bohr)')
ax3.set_ylabel('frequency')
ax3.grid()
ax3.set_title('Cumulative Psi and Psi^2 for {} Runs'.format(nReps))
ax3.legend(loc='best')

plt.show()
