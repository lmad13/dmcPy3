import numpy as np
import matplotlib.pyplot as plt
#DMC procedural programing style

#Define simulation parameters:
#Parameters below are simulation specific
#The length of time each timeStep is meant to simulate:
timeStep=10.0

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
nDimensions=1
mass=np.zeros(nDimensions*nAtoms) 
#mass=np.array([12.00000,15.995])/(6.02213670000e23*9.10938970000e-28)
#mass=np.array([10.00000,10.0000])/(6.02213670000e23*9.10938970000e-28)
#reducedMass=np.prod(mass)/np.sum(mass)
#reducedMass=10000.0
reducedMass=9114.44257
mass[:]=reducedMass
equilibriumPosition=5.0

print("reduced Mass", reducedMass)
print("atomic mass",2*reducedMass*(6.02213670000e23*9.10938970000e-28))
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
	
	#The bond length, x is:
	bondLength=np.linalg.norm(atomicPositions[:,0,:]-atomicPositions[:,1,:],axis=1)
	#x=bondLength-(1.1283*0.529177)
	
	x=bondLength-equilibriumPosition
	
	#The spring constant is:
	#omega=2169.8
	#omega=1000.0
	#k=(omega*2.0*np.pi*2.998*10**(10))**2*reducedMass #cm-1

	#convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
	#k=k*convfactor   #Eh/(bohr^2)
	
	k=1.0
	#print("k: ",k)
	#And the potential energy as a function of the spring constant and the bond lenght is:
	potentialEnergies=0.5*k*x**2

	return potentialEnergies
def HarmonicWavefunction(x,m):
	#omega=2169.8
	#k=(omega*2.0*np.pi*2.998*10**(10))**2*reducedMass #cm-1
	#convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
	#k=k*convfactor   #Eh/(bohr^2)
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

	for atom in range(nAtoms):
		for d in range(nDimensions):
			movement[:,atom,d]=np.random.normal(centerOfGaussian,widthOfGaussian[atom],size=(currentPopulation))


	newAtomicPositions=atomicPositions+movement


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
	
		#Diffuse the Walkers
		# if iStep in [100,500,1000]:	
		# 	oldBondLength=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		# 	oldEnergy=CalcPotentialEnergy(positions)
		# 	plt.scatter(oldBondLength,oldEnergy,marker='v')


		positions=Diffuse(positions)
		# if iStep in [100,500,1000]:	
		# 	newBondLength=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		# 	newEnergy=CalcPotentialEnergy(positions)
		# 	dE=newEnergy-oldEnergy
		# 	dBL=newBondLength-oldBondLength

		# 	plt.quiver(oldBondLength,oldEnergy,dBL,dE*10)

			
		# 	plt.plot([newBondLength,oldBondLength],[newEnergy,oldEnergy])
		# 	plt.scatter(newBondLength,newEnergy,marker='v')

		#Determine the potential Energy
		potentialEnergy=CalcPotentialEnergy(positions)
		bondLength=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		
		#plt.scatter(bondLength,potentialEnergy)
		#plt.title(str("LRM population:"+ str(currentPopulationSize)+ "  Reference Energy: "+str(referenceEnergy) ) )
		#plt.show()
	
		#Draw a random number from a normal distribution limited to be between 0 and 1.
		randomNumber=np.random.uniform(size=currentPopulationSize)
	
		#Determine the probability that the walkers will be deleted. 
		probabilityDelete=np.exp(-(potentialEnergy-referenceEnergy)*timeStep)
		#probabilityDelete=1.0-np.exp(-(potentialEnergy-referenceEnergy)*timeStep)

		#Compare the probability to the random number.  Those with a random number higher than the probability of deletion are NOT delete. 
		#Those walker survive.
		maskSurvive = (randomNumber<=probabilityDelete)
		maskDelete= (randomNumber>probabilityDelete)
	
		#Determine the probability of that the walkers are replicated.
		probabilityReplicate=np.exp(-(potentialEnergy-referenceEnergy)*timeStep)-1.0
	
		#Compare the probability to the random number.  Those with a random number lower than the probability of replication will be duplicated. 
		maskReplicate = (randomNumber<=probabilityReplicate)
		
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
	


		
	print("Average Bond Length" ,np.average(bondLengthList[int(nTimeSteps/2):])/0.529177)
	print("Std Bond Length" ,np.std(bondLengthList[int(nTimeSteps/2):])/0.529177)

	print("Zero Point Energy:",np.average(zeroPointEnergy[int(nTimeSteps/2):])*au2wn)
	#Do something with the data
	print("Theoretical harmonic:",2169.8/(2)," cm^-1", 2169.8/(au2wn*2),"Hartrees")
	print("Theoretical harmonic:",np.sqrt(1.0/reducedMass)/2*(au2wn)," cm^-1", np.sqrt(1.0/reducedMass)/2,"Hartrees")

	#plt.plot(np.arange(nTimeSteps),np.array(zeroPointEnergy)*au2wn)
	#plt.show()
	#plt.savefig("ZeroPointEnergy.png")
	


	return(positions,positionList,zeroPointEnergy)

#plt.plot(np.arange(nTimeSteps),populationList)
#plt.show()
#plt.savefig("population.png")

#Initialize Walkers
#Create nWalkers worth of replications of the molecular system.

#Each walker has an associated "position" represented by the x,y,z Cartesian
#coordinate of each atom in the molecule.

positions=np.zeros((nWalkers,nAtoms,nDimensions))
positions[:,0,0]=5.0

ZPE=[]
positions,positionList,refEnergy=Propagate(positions,10000)
for irep in range(5):
	positions,positionList,refEnergy=Propagate(positions,nTimeSteps)
	ZPE.append(np.average(refEnergy[int(nTimeSteps/2):]))
ZPE=np.array(ZPE)
print(ZPE*au2wn)
print(np.average(ZPE)*au2wn, np.average(ZPE))

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






