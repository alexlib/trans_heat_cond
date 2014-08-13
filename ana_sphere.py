# Analytical solution of 1D transient heat conduction in a sphere
# convection at surface & no heat generation
# see Incropera2011 7th Edition book, pg.303

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries and packages
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as py

#--- parameters from Papadikis2010a Table 1
rhow = 700      # density of biomass, 700 kg/m3
d = 0.035e-2    # particle diameter for 350um, m
cpw = 1500      # specific heat capacity biomass, J/(kg K)
kw = 0.105      # thermal conductivity biomass, W/(m K)
Ti = 300        # uniform initial temp of sphere, K
Tinf = 773      # surrounding fluid or gas temp, K
tmax = 0.8      # max time, s

#--- theta Calculations for particle diameter 350um
h = 375                 # heat transfer coefficent, W/(m2 K)
ro = (d/2)              # radius of sphere, m

alpha = kw/(rhow*cpw)   # thermal diffusivity biomass, m2/s

Bi = (h*ro)/kw          # Biot number, (-)

#--- function f(z) = 1 - z*cot(z) - Bi, note that cot = 1/tan

def funcZeta(z):
    f = 1 - z*(1 / np.tan(z)) - Bi
    return f

#--- find sign change and roots of function
step = 0.1
zi = 0
zf = 1250 + step  # add step to end on integer

f = []
root = []

for i in np.arange(zi, zf, step):
    f.append(funcZeta(i))

k = 0
s = 0

for i in range(0, len(f)):
    k = k + 1
    s = s + 0.1
    if f[k-1] < 0 and f[k] > 0:
        root.append(opt.brentq(funcZeta, s-0.1, s))

#--- calculate Fourier number
n = len(root)
t = np.arange(0, tmax+0.01, 0.01)
Fo = (alpha * t) / (ro**2)   # Fourier number, (-)

#--- for ro (outer surface)
r = ro

def funcCn(root):
    Cn = 4*(np.sin(root)-root*np.cos(root)) / (2*root-np.sin(2*root))
    return Cn

def funcDn(ro, r, root):
    Dn = (ro/(root*r)) * np.sin(root*(r/ro));
    return Dn

theta_o = funcCn(root[0])*np.exp(-root[0]**2 * Fo)*funcDn(ro,r,root[0])

for i in range(1, n):
    dTheta_o = funcCn(root[i])*np.exp(-root[i]**2 * Fo)*funcDn(ro,r,root[i])
    theta_o = theta_o + dTheta_o

# convert dimensionless theta to temperature K for ro (outer surface)
T_o = Tinf + theta_o*(Ti-Tinf)

#--- for r (center)
r = 1e-10

theta_r = funcCn(root[0])*np.exp(-root[0]**2 * Fo)*funcDn(ro,r,root[0])

for i in range(1, n):
    dTheta_r = funcCn(root[i])*np.exp(-root[i]**2 * Fo)*funcDn(ro,r,root[i])
    theta_r = theta_r + dTheta_r

# covert dimensionless theta to temperature K for r (center)
T_r = Tinf + theta_r*(Ti-Tinf)

#--- plot results

# configuration for plot
if Ti > Tinf:
    # cooling sphere, i.e. Ti=773K and Tinf=300K
    Th = Ti
    ylim =[Ti+20, Tinf-20]
else:
    # heating sphere, i.e. Ti=300K and Tinf=773K
    Th = Tinf
    ylimRange = [Ti-20, Tinf+20]

# plot
py.figure(1)
py.plot(t, T_o, '-r', label='surface')
py.plot(t, T_r, '--r', label='center')
py.title('Sphere (analytical)')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.ylim(ylimRange)
py.xlim([0, tmax])
py.axhline(y=Th, color='k', linestyle='--')
py.rcParams['xtick.major.pad'] = 6
py.rcParams['ytick.major.pad'] = 6
py.legend(loc='best', numpoints=1)
py.grid()
py.show()
