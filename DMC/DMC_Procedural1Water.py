import numpy as np
import matplotlib.pyplot as plt
import FlexibleSPCEPotentialComplete as PES # my file 
from scipy.optimize import curve_fit
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
nAtoms=3
nDimensions=3
mass=np.zeros(nDimensions*nAtoms) 
#mass=np.array([12.00000,15.995])/(6.02213670000e23*9.10938970000e-28) #CO
mass=np.array([15.995,1.00000,1.00000])/(6.02213670000e23*9.10938970000e-28) # (Oxyfen, hydrogen, hydrgen) ... but change mass matrix if we add more molecules 
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
	#For single water only
	return	PES.PotentialEnergySingleWater(atomicPositions) # might do PotentialEnergyTwoWaters


def HarmonicWavefunction(x,N,k,m,re):
	#this is the theoretical wavefunction equation.  A bit awkward for sure with the dependence on k that 
	#is defined in the potential energy function.  
	return N*np.exp(-((x-(re))**2)/2*np.sqrt(m*k)) # Gaussian function to fit histograms (how close is our histogram to the normal distribution)

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
	populationList=[] # how many walkers we have 


	for iStep in range(nTimeSteps):
	
		#Determine the Reference energy
		currentPopulationSize=positions.shape[0]
		averageEnergy=np.average(CalcPotentialEnergy(positions)) # zero point energy
		referenceEnergy=averageEnergy+ ((1.0-float(currentPopulationSize/nWalkers)) / (2.0*timeStep))
		# run with Python 3

		# diffuse the positions
		positions=Diffuse(positions) 
		
		#!Here is another example of the sort of plotting that I do to figure out if the diffusion,
		#potential energy functions are doing what I expect.  

		# print statements and plot statements
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
	
		#Determine the probability that the walkers will be DELETED. (~10)
		probabilityDelete=np.exp(-(potentialEnergy-referenceEnergy)*timeStep)

		#Compare the probability to the random number.  Those with a random number higher than the probability of deletion are NOT delete. 
		#Those walker survive.
		maskSurvive = (randomNumber<=probabilityDelete)
		maskDelete= (randomNumber>probabilityDelete)
	
		#Determine the probability of that the walkers are REPLICATED. (~10)
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
	print("Average, Median, StD of Pop List",np.average(populationList),np.median(populationList),np.std(populationList))

	return(positions,positionList,zeroPointEnergy)


#Initialize Walkers
#Create nWalkers worth of replications of the molecular system.

#Each walker has an associated "position" represented by the x,y,z Cartesian
#coordinate of each atom in the molecule.

positions=np.zeros((nWalkers,nAtoms,nDimensions)) # start initial positions at zero (but for more than 1 water, the atoms in another water are stacked up at another coordinates)
#start all of the walkers with the 1st atom (Oxygen) as position x=0, y=0, z=0 and 
# the second atom (hydrogen) at x=1, y=0, z=0 and 
# the third atom (hydrogen) at x=-1, y=-1, z=0
positions[:,1,0]=1.0
# positions[:, atom, coordinate] such that 1 in atom is first hydrogen and 0 in coordinate is x 
positions[:,2,0]=-1.0
positions[:,2,1]=-1.0
 

ZPE=[]
positionListCollection = []
positions,positionList,refEnergy=Propagate(positions,50) # equilibriate 
nRep=3 # repeat for 5 times of nTimeSteps (100) <-- can increase to 1000
for irep in range(nRep):
	positions,positionList,refEnergy=Propagate(positions,nTimeSteps)
	ZPE.append(np.average(refEnergy[int(nTimeSteps/2):]))
	positionListCollection.append(np.array(positionList))
	plt.ylim(.05,.07)
	plt.plot(np.arange(irep*nTimeSteps,irep*nTimeSteps+nTimeSteps),refEnergy,linewidth=.5)


ZPE=np.array(ZPE)

print("here's your ",nRep," energies: ",ZPE)
print("The average: ",np.average(ZPE))
print("The standard deviation: ",np.std(ZPE))

plt.hlines(np.average(ZPE),0,nRep*nTimeSteps)
plt.ylabel("Reference Energy (Hartree or a.u.)")
plt.xlabel("time step")
plt.show()

wavefunctionAsym=[]
print(positionListCollection[0])
for ensemble in positionListCollection:
	for positions in ensemble:
		# measure asymmetric stretch for water (if equal, value of stretch=0, if unequal, the value=either + or -)
		asymStretch=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)-np.linalg.norm(positions[:,0,:]-positions[:,2,:],axis=1)
		# histogram of how common the strecth value is (expect equal (0) to be most common (meaning the distribution is centered at 0))
		BLhistogram,bin_edges=np.histogram(asymStretch,bins=50,range=(-1.0,1.0))
		wavefunctionAsym.append(BLhistogram)
histogramXaxis=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.scatter(histogramXaxis,np.average(wavefunctionAsym,axis=0)/np.sum(np.average(wavefunctionAsym,axis=0)))
popt, pcov = curve_fit(HarmonicWavefunction, histogramXaxis, np.average(wavefunctionAsym,axis=0)/np.sum(np.average(wavefunctionAsym,axis=0)),
	bounds=(0, [1.0,500., 100., 0.05]))
print(popt)
print(1059.162 *(1.0/0.529177)**2 * (4.184/2625.5))

plt.plot(histogramXaxis, HarmonicWavefunction(histogramXaxis, *popt), 'r-',
	label='fit: N=%5.3f, a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.show()
wavefunctionBL=[]
for ensemble in positionListCollection:
	for positions in ensemble:
		BL1=np.linalg.norm(positions[:,0,:]-positions[:,1,:],axis=1)
		# for just one bond lenghth (distribution should be centered around the bond length)
		BLhistogram,bin_edges=np.histogram(BL1,bins=50,range=(1.00,3.0))
		wavefunctionBL.append(BLhistogram)
histogramXaxis=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.scatter(histogramXaxis,np.average(wavefunctionBL,axis=0)/np.sum(np.average(wavefunctionBL,axis=0)))
popt, pcov = curve_fit(HarmonicWavefunction, histogramXaxis, np.average(wavefunctionBL,axis=0)/np.sum(np.average(wavefunctionBL,axis=0)),
	bounds=(0, [1.0,200., 20., 2.0]))
print(popt)
print(1059.162 *(1.0/0.529177)**2 * (4.184/2625.5))

plt.plot(histogramXaxis, HarmonicWavefunction(histogramXaxis, *popt), 'r-',
	label='fit: N=%5.3f, a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.show()
wavefunctionAngle=[]
for ensemble in positionListCollection:
	# bond angle 
	aHOH=[]
	for walkerPos in ensemble:
		vecOH_1=walkerPos[0]-walkerPos[1]
		vecOH_2=walkerPos[2]-walkerPos[0]
		cosAngle=np.dot(vecOH_1,vecOH_2)/(np.linalg.norm(vecOH_1)*np.linalg.norm(vecOH_2))
		aHOH.append(np.arccos(cosAngle))

	aHOH=np.array(aHOH)
	aHOHeq= 112.0 * np.pi/180.0

	# expect the histogram to be centered about the ^ hard coded aHOHeq
	AngleHistogram,bin_edges=np.histogram(aHOH,bins=100,range=(0,np.pi))
	wavefunctionAngle.append(AngleHistogram)
histogramXaxis=(bin_edges[:-1]+bin_edges[1:])/2.0
plt.scatter(histogramXaxis,np.average(wavefunctionAngle,axis=0)/np.sum(np.average(wavefunctionAngle,axis=0)))
popt, pcov = curve_fit(HarmonicWavefunction, histogramXaxis, np.average(wavefunctionAngle,axis=0)/np.sum(np.average(wavefunctionAngle,axis=0)),
	bounds=(0, [1.0,500., 100., 2.0]))
print(popt)
print(1059.162 *(1.0/0.529177)**2 * (4.184/2625.5))

plt.plot(histogramXaxis, HarmonicWavefunction(histogramXaxis, *popt), 'r-',
	label='fit: N=%5.3f, a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

plt.show()



