import numpy as np
import matplotlib.pyplot as plt
#graph psi, psi squared, potential energy surface of Morse fit using hardcoded in values

#potential energy surface
scfdone='5.02685721e-03 3.10285933e-03 1.80171111e-03 9.53720631e-04 4.32416217e-04 1.43875692e-04 1.86496856e-05 5.65344352e-06 6.75517276e-05 1.77275682e-04 1.77275682e-04 3.15397403e-04 4.68153983e-04 6.25962970e-04 7.82309354e-04 9.32913123e-04 1.07510848e-03 1.20738251e-03 1.32903373e-03 1.43992072e-03 1.54027814e-03 1.63058310e-03 1.71145913e-03 1.78360790e-03 1.84776151e-03 1.90464998e-03'

scfdone2 = str.split(scfdone)
scfdone3=[]
for i in scfdone2:
    scfdone3.append(float(i))
scfdone = np.array(scfdone3)

#potential energy surface xdata
pxdata = '2.96 3.06 3.16 3.26 3.36 3.46 3.56 3.76 3.86 3.96 4.06 4.16 4.26 4.36'

pxdata2 = str.split(pxdata)
pxdata3 = []
for f in pxdata2:
    pxdata3.append(float(f))
pxdata3 = np.array(pxdata3)


#ydata
ydata = '5.080923e-03 3.048065e-03 1.743530e-03 9.284800e-04 4.387520e-04 1.643530e-04 3.343000e-05 0.000000e+00 3.592500e-05 1.222490e-04 1.222490e-04 2.469230e-04 3.990110e-04 5.686370e-04 7.457970e-04 9.228480e-04 1.090627e-03 1.244252e-03 1.379021e-03 1.495127e-03 1.590694e-03 1.668095e-03 1.728700e-03 1.775436e-03 1.811979e-03 1.839491e-03'

ydata2 = str.split(ydata)
ydata3 = []
for e in ydata2:
    ydata3.append(float(e))
ydata = np.array(ydata3)

#b7
b7 = '2.96 3.06 3.16 3.26 3.36 3.46 3.56 3.66 3.76 3.86 3.86 3.96 4.06 4.16 4.26 4.36 4.46 4.56 4.66 4.76 4.86 4.96 5.06 5.16 5.26 5.36'

b72 = str.split(b7)
b73=[]
for a in b72:
    b73.append(float(a))
b7 = np.array(b73)

#psi
psi = '4.99810086e-05 7.49925123e-05 1.34955030e-04 2.09969535e-04 2.74958043e-04 4.49896087e-04 6.69797648e-04 7.84780668e-04 1.12470524e-03 1.63456832e-03 2.15945493e-03 2.77433751e-03 3.38910122e-03 4.78374748e-03 5.71862809e-03 6.82324135e-03 8.62776320e-03 1.05873966e-02 1.28267002e-02 1.45463914e-02 1.65106618e-02 1.96048545e-02 2.11394943e-02 2.36041391e-02 2.63581084e-02 2.80279276e-02 3.02321735e-02 3.26066466e-02 3.33415031e-02 3.59361326e-02 3.61057197e-02 3.67207914e-02 3.72055835e-02 3.76303164e-02 3.80104311e-02 3.66258166e-02 3.53962719e-02 3.53359201e-02 3.36413741e-02 3.08071357e-02 2.97024929e-02 2.78879041e-02 2.63131708e-02 2.35740191e-02 2.19993684e-02 2.09295922e-02 1.85502952e-02 1.68554011e-02 1.50309995e-02 1.35014653e-02 1.14469838e-02 1.05773381e-02 9.06773827e-03 7.98291016e-03 6.85820585e-03 6.01344272e-03 4.98384439e-03 4.40885736e-03 3.73903974e-03 3.07421458e-03 2.63441350e-03 2.21441895e-03 1.83957583e-03 1.34967576e-03 1.21971873e-03 1.05974770e-03 8.49777175e-04 6.74855622e-04 5.74850103e-04 5.19871606e-04 3.79898073e-04 3.39920564e-04 2.34921058e-04 2.54918056e-04 2.14945549e-04 1.34954029e-04 7.49925085e-05 8.99745198e-05 5.49790143e-05 3.99930068e-05 2.99970034e-05 3.99910080e-05 3.99915070e-05 5.00050005e-06 9.99850225e-06 1.99910046e-05 0.00000000e+00 0.00000000e+00 1.49960045e-05 5.00000000e-06 0.00000000e+00 0.00000000e+00 5.00050005e-06 9.99900170e-06 4.99900020e-06 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 4.99750125e-06'

psi2 = str.split(psi)

psi3 = []
for b in psi2:
    psi3.append(float(b))

psi = np.array(psi3)

#psi squared
psisqrt = '1.49329569e-07 3.92188267e-07 9.10030901e-07 1.67785669e-06 3.11637673e-06 6.96978912e-06 1.57667891e-05 1.88659113e-05 4.18529640e-05 9.38316279e-05 1.59237983e-04 2.65450604e-04 4.12082422e-04 8.09952745e-04 1.15855458e-03 1.65966678e-03 2.68393566e-03 4.06225361e-03 6.06564140e-03 7.79790133e-03 1.01141275e-02 1.43223797e-02 1.68115976e-02 2.09252749e-02 2.62211903e-02 2.96616844e-02 3.46399239e-02 4.04245131e-02 4.22891795e-02 4.92970818e-02 4.97372497e-02 5.14913736e-02 5.27616779e-02 5.38323919e-02 5.49360669e-02 5.10053827e-02 4.75708714e-02 4.72313816e-02 4.28792041e-02 3.57834036e-02 3.32027923e-02 2.92780948e-02 2.60225627e-02 2.07876819e-02 1.81373003e-02 1.63744811e-02 1.27602543e-02 1.05605309e-02 8.34191180e-03 6.72354289e-03 4.82447148e-03 4.10421317e-03 3.02275443e-03 2.34067677e-03 1.71039146e-03 1.32766158e-03 9.23463737e-04 6.93880183e-04 5.08497451e-04 3.40947829e-04 2.54327104e-04 1.79061921e-04 1.22217589e-04 6.75415515e-05 5.98274808e-05 4.16580792e-05 2.81837941e-05 1.97251503e-05 1.43624228e-05 1.11172342e-05 6.33365305e-06 5.16528858e-06 2.43607890e-06 2.94053195e-06 1.78755979e-06 8.77769541e-07 5.53872973e-07 3.72205720e-07 3.57267027e-07 1.68489513e-07 1.29626342e-07 1.85995207e-07 1.30629988e-07 1.86097753e-08 3.72486461e-08 7.48560584e-08 0.00000000e+00 0.00000000e+00 5.62750220e-08 1.86685933e-08 0.00000000e+00 0.00000000e+00 0.00000000e+00 3.73018055e-08 1.88979902e-08 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00'

psisqrt2 = str.split(psisqrt)
psisqrt3 = []
for c in psisqrt2:
    psisqrt3.append(float(c))
psisqrt = np.array(psisqrt3)

#bin
bin = '2.97972574 3.00041244 3.02109915 3.04178585 3.06247255 3.08315925 3.10384595 3.12453265 3.14521935 3.16590605 3.18659275 3.20727946 3.22796616 3.24865286 3.26933956 3.29002626 3.31071296 3.33139966 3.35208636 3.37277307 3.39345977 3.41414647 3.43483317 3.45551987 3.47620657 3.49689327 3.51757997 3.53826667 3.55895338 3.57964008 3.60032678 3.62101348 3.64170018 3.66238688 3.68307358 3.70376028 3.72444698 3.74513369 3.76582039 3.78650709 3.80719379 3.82788049 3.84856719 3.86925389 3.88994059 3.9106273  3.931314   3.9520007 3.9726874  3.9933741  4.0140608  4.0347475  4.0554342  4.0761209 4.09680761 4.11749431 4.13818101 4.15886771 4.17955441 4.20024111 4.22092781 4.24161451 4.26230121 4.28298792 4.30367462 4.32436132 4.34504802 4.36573472 4.38642142 4.40710812 4.42779482 4.44848153 4.46916823 4.48985493 4.51054163 4.53122833 4.55191503 4.57260173 4.59328843 4.61397513 4.63466184 4.65534854 4.67603524 4.69672194 4.71740864 4.73809534 4.75878204 4.77946874 4.80015544 4.82084215 4.84152885 4.86221555 4.88290225 4.90358895 4.92427565 4.94496235 4.96564905 4.98633576 5.00702246 5.02770916'

bin2 = str.split(bin)
bin3 = []
for d in bin2:
    bin3.append(float(d))
bin = np.array(bin3)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)

ax1.plot(bin, psi, label = 'Psi')
ax1.plot(bin, psisqrt, label = 'Psi Squared')
ax1.set_xlim(np.amin(bin),np.amax(bin))
ax1.set_xlabel('Oxygen-Carbon Bond Length (Angstrom)',fontsize=17)
ax1.set_ylabel('Frequency',fontsize=17)
ax1.tick_params(labelsize=13)
ax1.grid()
ax1.legend(fontsize=15)


ax2.plot(b7, scfdone*220000.00000000003,'--', label = 'Morse Fit')
ax2.plot(b7, ydata*220000.00000000003, 'o', label = 'SCF Energy')
ax2.set_xlim(np.amin(bin),np.amax(bin))
ax2.axhline(y=0.00019319303655417793*220000.00000000003,c='red',label='Zero-point energy')
ax2.axhline(y=0.00028512262618178233*220000.00000000003,c='blue',label='Zero-point energy for Harmonic')
ax2.legend(fontsize=15)
ax2.grid(axis='x')
ax2.tick_params(labelsize=13)
ax2.set_ylabel('SCF Energy (cm^-1)',fontsize=17)
ax2.set_xlabel('Oxygen-Carbon Bond Length (Angstrom)',fontsize=17)
plt.show()

