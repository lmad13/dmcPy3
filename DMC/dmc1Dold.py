
"""
This is the dmc code I obtained from my GitHub library that contains all the old code. Because harmonicOscillatorH2_5
stopped working somehow, I tried to see if this fixes anything. I added the probability of crossing over to thise code.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


au2wn=219474.63

class Wavefunction:
    
    def __init__(self, nWalkers, potential, plotting=False, dtau= 10.0, omegaInput=0.00057161, molecule1='H2', molecule2='H2O', molecule3='H2O', recrossing=True):
        """
        Wavefunction constructor.
        Creates ...
        
        Parameters:
        ---------------
        nWalkers: int. Number of walkers
        potential: str. Potential energy function (harmonic, half harmonic right, half harmonic left, Morse, Quartic, half harmonic)
        plotting: bool. Whether not to plot (DEFAULT: False)
        omegaInput: flt. Omega number for k (spring constant) (DEFAULT: 0.00057161)
        dtau: flt. time step (DEFAULT: 10.0)
        molecule1: str. atom/molecule 1 (DEFAULT: H2)
        molecule2: str. atom/molecule 2 (DEFAULT: H20)
        molecule3: str. atom/molecule 3 (DEFAULT: H20)
        recrossing: bool. Choose to recross (DEFAULT: True)
        """

        # print statement for how many walkers and potential energy surface
        print("initialized"+str(nWalkers)+" coordinates for " +str(potential))

        # x coordinates
        self.xcoords=np.zeros(nWalkers)
        # set potential energy surface and print it out
        self.potential=self.setPotential(potential)
        # set time step and print it out
        self.dtau=self.set_dtau(dtau)
        # set reduced mass and print it out
        self.mass=self.set_mass(molecule1, molecule2, molecule3)
        # set plot feature and print it out
        self.plotting=self.setPlotting(plotting)
        # set recrossing feature and print it out
        self.recrossing=self.setRecrossing(recrossing)
        
        self.D=0.5
        # Omega number for k (spring constant)
        self.omega=omegaInput
        self.sigma_dx=(2.0*self.D*self.dtau/self.mass)**0.5
        print('sigma_dx'+str(self.sigma_dx))
        self.mu_dx=0.0
        self.alpha=0.500/self.dtau


    # get and set functions
    def setPotential(self,potential):
        """
        Set potential energy from harmonic, half harmonic right, half harmonic left, Morse, Quartic, half harmonic
        and print what it is

        Parameters:
        ---------------
        potential: str. Potential energy function (harmonic, half harmonic right, half harmonic left, Morse, Quartic, half harmonic)
        """
        print('potential surface is '+ str(potential))
        return potential


    def setPlotting(self, plotting):
        """
        Set whether you want to plot or not
        and print what it is

        Parameters:
        ---------------
        plotting: bool. whether you want to plot or not
        """
        print('plotting is '+ str(plotting))
        return plotting


    def setRecrossing(self, recrossing):
        """
        Set whether you want to recross or not
        and print what it is

        Parameters:
        ---------------
        recrossing: bool. whether you want to recross or not
        """
        print('recrossing is '+ str(recrossing))
        return recrossing
        

    def setX(self, x):
        """
        Set xcoords instance variable to x

        Paramters:
        ------------------
        x: int. number of walkers
        """
        self.xcoords=x


    def set_dtau(self, dtau):
        """
        Set dtau instance variable as dtau

        Paramters:
        ------------------
        dtau: flt. time step
        """
        print('set dtau to be '+ str(dtau))
        return dtau


    def get_dtau(self):
        """
        Get dtau instance variable
        """
        return self.dtau

    def return_mass(self, molecule):
        """
        Obtain molecule and return reduced mass
        Available: 'H2', 'Li2' , 'H2O', 'HO', 'CH4'    

        Paramters:
        -------------
        molecule: str. Atomic 
        """
        #why? for reduced mass in amu
        conversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839    in Atomic Units!!
        
        # mass of atoms used
        massH=1.00782503223
        massLi=7.0
        massO=15.994915
        massC=12.000000

        # initialize massAtom
        massAtom = 0

        # set massAtom depending on the molecule
        if molecule=='H2':
            massAtom=massH*2
        elif molecule=='Li2':
            massAtom=massLi*2
        elif molecule=='H2O':
            massAtom=massO+massH*2
        elif molecule=='HO':
            massAtom=massO+massH
        elif molecule=='CH4':
            massAtom=massH*4+massC

        return massAtom


    def set_mass(self, molecule1, molecule2, molecule3):
        """
        Get reduced mass between different molecules
        Available options: 'H2', 'Li2' , 'H2O', 'HO', 'CH4' 

        Paramters:
        --------------------
        molecule1: str. atom/molecule 1 (DEFAULT: H2)
        molecule2: str. atom/molecule 2 (DEFAULT: H20)
        molecule3: str. atom/molecule 3 (DEFAULT: H20)
        """

        #why? for reduced mass in amu
        conversionFactor=1.000000000000000000/6.02213670000e23/9.10938970000e-28#1822.88839    in Atomic Units!!
        
        # mass of atoms used
        massH=1.00782503223
        massLi=7.0
        massO=15.994915
        massC=12.000000
        
        # obtain reduced mass of each molecule
        massAtom1 = self.return_mass(molecule1)
        massAtom2 = self.return_mass(molecule2)
        massAtom3 = self.return_mass(molecule3)

        # reduced mass with conversino factor
        reducedMass = (massAtom1*massAtom2*massAtom3)/(massAtom1+massAtom2+massAtom3)*conversionFactor
        # print statement for reduced mass
        print('reduced mass is '+str(reducedMass))
        
        return reducedMass


    def getTheoreticalOmega0(self):
        """

        """
        if self.potential=='harmonic':
            return self.omega/(2.0*au2wn)
        
        elif 'half harmonic' in self.potential:
            return 3.0*self.omega/(2.0*au2wn)
        
    def getAlpha(self):
        k=(self.omega*2.0*np.pi*3.0*10**(10))**2*self.mass #cm-1                                            
        convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2   
        k=k*convfactor

        alpha=np.sqrt(k*self.mass)
        return alpha

    def crossProb(self, originalx, newx):
        insideNumerator = -2*self.mass*originalx*newx
        output = np.exp(insideNumerator/self.dtau)  
        return output      


    def V(self,x):
        #input in bohr, output in a.u.
        #k=(self.omega*2.0*np.pi*3.0*10**(10))**2*self.mass #cm-1
        #convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
        #k=k*convfactor
        k = (self.omega**2)*self.mass
        #input for Morse
        De=0.0023164143974554563
        a=0.7228166340696739
        re=6.8463039861977375
        
        qa=0.09223635938249808
        qb=-2.973531354477643*10**(-8)
        qc=0.0747096550705948
        qd=3.178170582263107*10**(-8)
        qe=0.00034241789678879364

        if self.potential=='harmonic':
            v=0.5*k*(x*x)
        if self.potential=='half harmonic right':
            #v=0.5*k*(x*x) if x>0 else 100000
            # inf=10000.0
            v=0.5*k*(x*x)
            mask=(x<0.0)
            v[mask]=100*0.5*k*(1.3915484473258797*1.3915484473258797)
            #v=[vb if c else vinf for (c,vb,vinf) in zip(x>0.0,0.5*k*(x*x),inf*np.ones(x.size))]
        if self.potential=='half harmonic left':
            #v=0.5*k*(x*x) if x>0 else 100000
            # inf=10000.0
            v=0.5*k*(x*x)
            mask=(x>0.0)
            v[mask]=100*0.5*k*(1.3915484473258797*1.3915484473258797)
        if self.potential=='Morse':
            v=De*((1-np.exp(-1.0*a*(x-re)))**2)
            
        if self.potential=='Quartic':
            v=qa*(x**4)+qb*(x**3)+qc*(x**2)+qd*(x)+qe
        elif self.potential=='half harmonic':
            v=0.5*k*(x*x)

        return v


    #function propagate
    def propagate(self,x,nSteps,setV_ref=False,ConstantV_ref=0,printCensus=False,nSize=0,plotWalkers=False,printFlag=False):

        descendants=np.zeros((x.size))  #not yet implemented
        whoYaFrom=np.arange(x.size)
        N_size_step=x.size
        if nSize==0:#if it is the default value of zero...it needs to be set.
            nSize=x.size
        if setV_ref:
            v_ref=ConstantV_ref
        else:
            v_ref=np.average(self.V(x))#+(self.alpha*(1-float(N_size_step)/float(nSize)))
        vRefList=[]
        population=[]
        if plotWalkers:
            plt.figure(2)

        for step in range(nSteps):
            dx=self.diffuse(x)
            originalx=x
            x=x+dx

            # crossProbability = np.exp(-2.0*(x-dx)*x*np.sqrt(self.mass*self.mass)/self.dtau) # this is the old probability of recrossing 
            crossProbability = self.crossProb(originalx, x) # has exponential function? this is the new probability of recrossing based on "Diffusion Monte Carlo in Internal Coordinates"
            # print(crossProbability)
            #plot of x, crossProbability
            #don't have it plt.plot and plt.show at the end of everything

            # calculate the pontential energy of new step
            v=self.V(x)



            if step%(500) and printFlag:
                print('step: '+ str(step))
                sys.stdout.flush()
            if self.plotting and step%(nSteps/10)==0:
                plt.scatter(x,v)
                plt.xlabel('Hydrogen-Water Bond Length (Angstrom)')
                plt.ylabel('Probability of Recross Death')
                xrange= np.linspace(np.min(x), np.max(x), num=5)
                plt.plot(xrange,np.zeros((5))+v_ref)
                plt.show()
            #Elimintation of a random selection of walkers in
            #the classically forbidden region
            N_r=np.random.random(N_size_step)
            #rewrite into single exponential 
            P_d=1-np.exp(-(v-v_ref)*self.dtau) # probability of death calculation (1-exponential function?)

            Diff=N_r-P_d
            mask_survive = (Diff>0)
            nDeaths=np.sum(np.array(Diff<0).astype(int))
            survivors=x[mask_survive]
            
            if printCensus: print('Census: Deaths:'+ str(nDeaths))

            #recrossing correction for pop near node, only true for half harmonic PES
            if self.recrossing:
                crossed=((x-dx)*x<0)
                # P_recrossDeath=np.exp(-2.0*(x-dx)*x*np.sqrt(self.mass*self.mass)/self.dtau) ##mass is reduced mass! # change this line to my function results
                P_recrossDeath = crossProbability # product of probability that the walker is too high + probability that the walker has crossed the nodal surface
                if self.plotting and np.sum(P_recrossDeath)!=0:
                    plt.scatter(x,P_recrossDeath)
                    plt.xlabel('Hydrogen-Water Bond Length (Angstrom)')
                    plt.ylabel('Probability of Recross Death')

                    plt.show()
                N_r_recross=np.random.random(N_size_step)
                Diff=N_r_recross-P_recrossDeath
                mask_survive_recross=(Diff>0.0000000000)
                tempRecrossCensus=np.sum(np.array(Diff<0).astype(int))
                mask_survive=np.logical_and(mask_survive, mask_survive_recross)
                survivors=x[mask_survive]            

            #Creation of a random selection of walkers in the classically allowed region
            P_exp_b=np.exp(-(v-v_ref)*self.dtau)-1.0
            P_exp_b[np.logical_not(mask_survive)]=0.000000 #BECAUSE THE DEAD CANNOT REPRODUCE!! NO ZOMBIE MOTHERS!! 
            if printCensus: print('\n P_exp_b'+ str(np.average(P_exp_b[(P_exp_b>0)]))+"  "+ str(np.std(P_exp_b[(P_exp_b>0)])))
            weight_P_b=P_exp_b.astype(int)
            P_b=P_exp_b-weight_P_b            #classically allowed region
            #P_b[np.logical_not(mask_survive)]=0.0
            Diff=N_r-P_b
            mask_b = (Diff<0)
            next_gen=x[mask_b] 
            new_pop_whoYaFrom=whoYaFrom[mask_b]
            nBirths=np.sum(np.array(Diff<0).astype(int))
            addBirthtot=0
            new_pop=next_gen
            #for the additional births
            for n,(particle,weight) in enumerate(zip(x,weight_P_b)):
                if weight>0: #i.e. the dead can't reproduce
                    if weight>10:
                        #this really shouldn't happen
                        print('weight is too big, resetting to 10. The weight is ',weight, ". The time step is ", step)
                        print(x[n],v[n],'<',v_ref, -(v[n]-v_ref))
                        weight=10
                    addBirthtot=addBirthtot+weight

                    temp=np.tile(particle,weight)
                    temp_whoYaFrom=np.tile(whoYaFrom[n],weight)
                    new_pop=np.concatenate((new_pop,temp))
                    new_pop_whoYaFrom=np.concatenate((new_pop_whoYaFrom,temp_whoYaFrom))

            if printCensus: print('. Births:'+ str(nBirths)+ ". Add' births: "+ str( addBirthtot))
            if plotWalkers:
                plt.scatter(x[mask_survive],v[mask_survive],c='black')
                plt.scatter(x[mask_b],v[mask_b],c='blue',s=(P_exp_b[mask_b]*100.0)**2)
                plt.scatter(x[(np.logical_not(mask_survive))],v[(np.logical_not(mask_survive))],c='red',s=(P_d[(np.logical_not(mask_survive))]*100.0)**2)
                plt.plot([-1.0,1.0],[v_ref,v_ref],c='magenta')
                plt.xlabel('Hydrogen-Water Bond Length (Angstrom)')
                plt.ylabel('Probability')
                #plt.quiver(x,v,dx,np.ones(N_size_step))
                plt.show()

            #readjust V_ref
            next_gen=new_pop
            next_gen_whoYaFrom=new_pop_whoYaFrom
            #collect survivors and next generation
            new_population=np.concatenate((survivors,next_gen))
            N_size_step=new_population.size
            v_average=np.average(self.V(new_population))
            if not setV_ref:
                v_ref=v_average+(self.alpha*(1.00-float(N_size_step)/float(nSize)))
            if printCensus: print('('+ str(N_size_step)+' / '+ str(nSize)+') v_ref '+ str(v_ref)+ ' = ' + str(v_average)+ ' + '+ str(self.alpha*(1-float(N_size_step)/float(nSize))))

            if v_ref<0 and step>5:
                print('this is problematic.  NSize is probably too small')
                print(' step:'+ str(step)+"  "+ str(v_ref)+ ' : '+ str(float(N_size_step)/float(nSize))+ ' = '+ str(float(N_size_step))+'/'+str(float(nSize)))

            vRefList.append(v_ref)
            population.append(N_size_step)
            whoYaFrom=np.concatenate((whoYaFrom[mask_survive],next_gen_whoYaFrom))
            x=new_population
            #print x.shape, whoYaFrom.shape
        #print 'Average from ',nSteps/2,' steps until the end is ',np.average(vRefList[nSteps/2:])*au2wn,
        #print 'final number of ancestors', N_size_step
        for anc in whoYaFrom:
            descendants[anc]=descendants[anc]+1

        return vRefList, population, x, descendants


    def diffuse(self,x):
        N_size_step=x.shape[0]
        dx=np.random.normal(self.mu_dx, self.sigma_dx, N_size_step)
        return dx
