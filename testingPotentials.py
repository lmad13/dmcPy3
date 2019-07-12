import dmc1D
import numpy as np
import matplotlib.pyplot as plt

def Morse(x, De=0.000628, Re=6.01, Alpha=0.777):
    V=De*(1.0-np.exp(-Alpha*(x)))**2
    return V

def Harmonic(x):
    omega=2000.0

    k=(omega*2.0*np.pi*3.0*10**(10))**2 #cm-1
    convfactor=9.10938291e-31*(1.0/4.35974417e-18)*(5.2917721092e-11)**2     #kg/amu Eh/J  m**2/Bh**2
    k=k*convfactor
    V=0.5*k*x*x
    return V


x=np.arange(-1.5, 1.50,.01)
morseV=Morse(x)
harmonicV=Harmonic(x)
massHH=dmc1D.wavefunction.set_mass(1,'H2')
print(massHH)
massH2H2O=dmc1D.wavefunction.set_mass(1,'H2H2O')
print(massH2H2O)

plt.plot(x*massHH,harmonicV)
plt.plot(x*massH2H2O,morseV)
plt.show()
