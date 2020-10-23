import numpy as np
import matplotlib.pyplot as plt
#DMC procedural programing style

#Define simulation parameters:
#Parameters below are simulation specific
#The length of time each timeStep is meant to simulate:
timeStep=1.0

#The number of walkers in the simulation
nWalkers=1000

#The number of timesteps in the simulation 
nTimeSteps=100

au2wn=219474.63

#The List of data that will be collected during the simulation.
#If you choose you could also write this data to disk, but I usually just gather it in RAM as the simulation runs.

#Parameters below are molecule specific.
#Each walker represents the molecule.  The Test molecule is H2.  It has 2 atoms
nAtoms=2
nDimensions=3
mass=np.zeros(nDimensions*nAtoms) 
#mass=np.array([12.00000,15.995])/(6.02213670000e23*9.10938970000e-28) #CO
mass=np.array([10.00000,10.0000])/(6.02213670000e23*9.10938970000e-28) #Fake Molecule ~B2
reducedMass=np.prod(mass)/np.sum(mass)

#reducedMass=9114.44257  #Setting this 
#mass[:]=reducedMass 	#For 1D system, mass is reduced Mass, for 3D+ mass is mass or atom
equilibriumPosition=5.0 #equilibrium bond length

print("reduced Mass", reducedMass)
print("atomic mass",2*reducedMass*(6.02213670000e23*9.10938970000e-28))

#printing/plotting controls
plotWalkers=False
printCensus=False

def CalcPotentialEnergy(atomicPositions):
	#Output the calculated potential energy in atomic energy units (Hartrees) 
	#Input:The positional cartesian coordinates or each atom of each walker in units of 
	#atomic length units (Bohr). Dimensionality of atomicPositions=[nWalkers,nAtoms, nDimensions]

	# For testing purposes I suggest you start with the potential known as the 
	#Harominc Oscillator for H2 Bond stretch:
	#V(x)= 1/2 k x^2
	#Where k is the spring constant for the vibrating bond and x is the bond length.
	#The spring constant for the H2 bond is 0.359540726 Atomic Unit of Energy/ 
	#(Atomic Unit of Length)^2)
	#This is a good function to start with for testing purposes because the 
	#reference energy for the harmonic oscillator is known to be:
	# E_ref = 1/2 * hbar * sqrt(k/m) 

	#The three lines below are specific to H2 potential energy function. 
	
	#The bond length is the distance between the two atoms:
	bondLength=np.linalg.norm(atomicPositions[:,0,:]-atomicPositions[:,1,:],axis=1)

	#variable x is difference between bond length and equilibrium bond length
	x=bondLength-equilibriumPosition
	
	#The spring constant is:
	k=1.0	
	#print("k: ",k)
	
	#The harmonic frequency (not actually used in the eqns, but useful for Chemists/Physcists to compare to experiment)

	omega=np.sqrt(k/reducedMass)/(2.0*np.pi*2.998*10**(10)) #cm-1
	
	#if you know omega (from experiment) you can determing k using:	
	#k=(omega*2.0*np.pi*2.998*10**(10))**2*reducedMass 
	#convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
	#k=k*convfactor   #Eh/(bohr^2)
	
	
	#And the potential energy as a function of the spring constant and the bond lenght is:
	potentialEnergies=0.5*k*x**2

	return potentialEnergies


def HarmonicWavefunction(x,m):
	#this is the theoretical wavefunction equation.  A bit awkward for sure with the dependence on k that 
	#is defined in the potential energy function.  

	k=1.0
	return np.exp(-((x-(equilibriumPosition))**2)/2*np.sqrt(reducedMass*k))

def Diffuse(atomicPositions):
	#The heart of the algorithm!  The atomic positions are randomly displaced.  The magnitude 
	#of the displacement of each atom is related to the mass of the atoms. The displacement is determined 
	#from a normal distribution (a.k.a. a Gaussian Function) of random numbers. The distribution is centered at zero 
	#and the width of the distribution is Sqrt(timeStep/mass)
	#Because the width of the distribution depends on the mass of the atom, we must iterate through the atoms as 
	#np.random.normal doesn't have an elegant way to change the width of the distribution otherwise.

	currentPopulation=atomicPositions.shape[0]
	movement=np.zeros(atomicPositions.shape)
	
	DiffusionConstant=0.50000
	centerOfGaussian=0.000000
	widthOfGaussian=np.sqrt(2.00000 * DiffusionConstant* timeStep/mass)

	for atom in range(nAtoms):  #for the 1D case, only 1 atom can move in 1 dimension!  #this would be different for 
							#the multiD case!
		for d in range(nDimensions):

			movement[:,atom,d]=np.random.normal(centerOfGaussian,widthOfGaussian[atom],size=(currentPopulation))


	newAtomicPositions=atomicPositions+movement

	#!Here's some example scratch code with how I figure out if my difffusion code is doing what I expect. 

	#newBondLength=np.linalg.norm(newAtomicPositions[:,1,:]-newAtomicPositions[:,0,:],axis=1)
	#oldBondLength=np.linalg.norm(atomicPositions[:,1,:]-atomicPositions[:,0,:],axis=1)
	#print("average bond length change: ",np.average(newBondLength-oldBondLength))
	#print(" std of bond length change: ",np.std(newBondLength-oldBondLength))
	# # print(bondLength)
	# BLhistogram,bin_edges=np.histogram(newBondLength-oldBondLength,bins=50)
	# histogramXaxis=(bin_edges[:-1]+bin_edges[1:])/2.0
	# plt.plot(histogramXaxis,BLhistogram)
	# plt.show()

	return newAtomicPositions


def Propagate(positions,nTimeSteps):
	zeroPointEnergy=[]
	positionList=[]
	bondLengthList=[]
	populationList=[]


	for iStep in range(nTimeSteps):
	
		#Determine the Reference energy
		currentPopulationSize=positions.shape[0]
		averageEnergy=np.average(CalcPotentialEnergy(positions)) 
		referenceEnergy=averageEnergy+ ((1.0-(currentPopulationSize/nWalkers)) / (2.0*timeStep))

		positions=Diffuse(positions)
		
		#!Here is another example of the sort of plotting that I do to figure out if the diffusion,
		#potential energy functions are doing what I expect.  

		# if iStep in [100,500,1000]:	
		# 	newBondLength=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		# 	newEnergy=CalcPotentialEnergy(positions)
		# 	dE=newEnergy-oldEnergy
		# 	dBL=newBondLength-oldBondLength
		#!I would only use one of these at a time
		# 	plt.quiver(oldBondLength,oldEnergy,dBL,dE*10)	
		# 	plt.plot([newBondLength,oldBondLength],[newEnergy,oldEnergy])
		# 	plt.scatter(newBondLength,newEnergy,marker='v')

		#Determine the potential Energy
		potentialEnergy=CalcPotentialEnergy(positions)
		bondLength=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		
		#!More sample plotting for what I look at when debugging.

		#plt.scatter(bondLength,potentialEnergy)
		#plt.title(str("LRM population:"+ str(currentPopulationSize)+ "  Reference Energy: "+str(referenceEnergy) ) )
		#plt.show()
	

		#Draw a random number from a normal distribution limited to be between 0 and 1.
		randomNumber=np.random.uniform(size=currentPopulationSize)
	
		#Determine the probability that the walkers will be deleted. 
		probabilityDelete=np.exp(-(potentialEnergy-referenceEnergy)*timeStep)

		#Compare the probability to the random number.  Those with a random number higher than the probability of deletion are NOT delete. 
		#Those walker survive.
		maskSurvive = (randomNumber<=probabilityDelete)
		maskDelete= (randomNumber>probabilityDelete)
	
		#Determine the probability of that the walkers are replicated.
		probabilityReplicate=np.exp(-(potentialEnergy-referenceEnergy)*timeStep)-1.0
	
		#Compare the probability to the random number.  Those with a random number lower than the probability of replication will be duplicated. 
		maskReplicate = (randomNumber<=probabilityReplicate)
		

		#! Here's even more plotting and printing examples of the data that I look at while debugging!
		#! I supress printing and decrease verbosity with the printCensus=False and plotWalker=False
		if printCensus and iStep%100==0:
			print("Census: step",iStep,"\n  Death: ", np.sum(maskDelete))
			print("  Briths:", np.sum(maskReplicate))
			print('('+ str(currentPopulationSize)+' / '+ str(nWalkers)+') v_ref '+ str(referenceEnergy)+ ' = ' + str(averageEnergy)+ ' + '+ str( (1.0-(currentPopulationSize/nWalkers)) / (2.0*timeStep)))

		if plotWalkers and iStep in [100,500,1000]:
			plt.title("Delta E vs probability")
			plt.xlabel("Delta E")
			plt.ylabel("probability")
			plt.scatter(potentialEnergy-referenceEnergy,probabilityDelete,color='red')
			plt.scatter(potentialEnergy-referenceEnergy,probabilityReplicate,color='blue')
			plt.show()

			plt.title("Bond length vs Probabilities")
			plt.xlabel("Bond length")
			plt.ylabel("probability")
			plt.scatter(bondLength,probabilityDelete,color='red')
			plt.scatter(bondLength,probabilityReplicate, color='blue')
			plt.show()

			plt.title("Bond Length vs Energy")
			plt.xlabel("Bond length")
			plt.ylabel("E")
			plt.scatter(bondLength[maskSurvive],potentialEnergy[(maskSurvive)],facecolors='none', edgecolors='black')
			plt.scatter(bondLength[maskDelete],potentialEnergy[(maskDelete)],facecolors='none',edgecolor='red',s=(probabilityDelete[maskDelete]*100.0)**2)
			plt.scatter(bondLength[maskReplicate],potentialEnergy[maskReplicate],facecolors='none',edgecolor='blue',s=(probabilityReplicate[maskReplicate]*100.0)**2)
			plt.hlines(referenceEnergy,np.min(bondLength),np.max(bondLength))
			plt.show()

		#Use the mask arrays to select the walkers that survive and the walkers that are replicated.  
		#Reassign positions of the walkers for the next round of the simulation.
		if np.sum(maskReplicate)>0:
			positions=np.concatenate((positions[maskSurvive],positions[maskReplicate]),axis=0)
		else:
			positions=positions[maskSurvive]
	
		#Collect some data for post simulation analysis (as applicable to your simulation question).
	
		zeroPointEnergy.append(referenceEnergy)
		positionList.append(positions)
		populationList.append(currentPopulationSize)
	
		bondLengthList.append(np.average(np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)))
	
	#!Sometimes it is handy to print the population as a function of time step. Here's how I do that.
	#plt.plot(np.arange(nTimeSteps),populationList)
	#plt.show()
	#plt.savefig("population.png")
	#plt.plot(np.arange(nTimeSteps),zeroPointEnergy)
	#plt.show()
	#plt.savefig("runEnergy.png")

	#Do something with the data, like print it out averages and standard deviations	
	print("Average Bond Length" ,np.average(bondLengthList[int(nTimeSteps/2):]))
	print("Std Bond Length" ,np.std(bondLengthList[int(nTimeSteps/2):]))
	print("Zero Point Energy:",np.average(zeroPointEnergy[int(nTimeSteps/2):]),np.std(zeroPointEnergy[int(nTimeSteps/2):]))
	print("Theoretical harmonic:", np.sqrt(1.0/reducedMass)/2,"Hartrees")

	return(positions,positionList,zeroPointEnergy)


#Initialize Walkers
#Create nWalkers worth of replications of the molecular system.

#Each walker has an associated "position" represented by the x,y,z Cartesian
#coordinate of each atom in the molecule.

positions=np.zeros((nWalkers,nAtoms,nDimensions))
#start all of the walkers with the 1st atom as position x=5, y=0, 0=0 and the second atom at x=0, y=0, z=0
positions[:,0,0]=5.0  

ZPE=[]
positions,positionList,refEnergy=Propagate(positions,1000)
nRep=5
for irep in range(nRep):
	positions,positionList,refEnergy=Propagate(positions,nTimeSteps)
	ZPE.append(np.average(refEnergy[int(nTimeSteps/2):]))
	plt.ylim(.003,.007)
	plt.plot(np.arange(irep*nTimeSteps,irep*nTimeSteps+nTimeSteps),refEnergy,linewidth=.5)


ZPE=np.array(ZPE)

print("here's your ",nRep," energies: ",ZPE)
print("The average: ",np.average(ZPE))
print("The standard deviation: ",np.std(ZPE))

plt.hlines(np.average(ZPE),0,nRep*nTimeSteps)
plt.ylabel("Reference Energy (Hartree or a.u.)")
plt.xlabel("time step")
plt.show()

wavefunction=[]

for ensemble in positionList[int(nTimeSteps/2):]:
	bondLength=np.linalg.norm(ensemble[:,0,:]-ensemble[:,1,:],axis=1)
	BLhistogram,bin_edges=np.histogram(bondLength,bins=50,range=(4.5,5.5))
	wavefunction.append(BLhistogram)
histogramXaxis=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.scatter(histogramXaxis,np.average(wavefunction,axis=0)/np.sum(np.average(wavefunction,axis=0)))
plt.plot(histogramXaxis,HarmonicWavefunction(histogramXaxis,reducedMass)/np.sum(HarmonicWavefunction(histogramXaxis,reducedMass)))

plt.show()

#Propagate the Walkers!






