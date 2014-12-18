"""
Analytical solution of 1D transient heat conduction in a solid sphere, cylinder,
or slab with convection at the surface and no heat generation.

References: 
1) Recktenwald 2006
2) Bergman, Lavine, Incropera, Dewitt 2011 from Ch. 5, pg.299-304
3) Papadikis 2010a
"""

# Modules
#------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as py
from funcTheta import theta

py.close('all')

# Parameters from Papadikis 2010a Table 1
#------------------------------------------------------------------------------

rhow = 700      # density of biomass, 700 kg/m3
d = 0.035e-2    # particle diameter for 350 um size, m
cpw = 1500      # specific heat capacity biomass, J/(kg K)
kw = 0.105      # thermal conductivity biomass, W/(m K)
Ti = 300        # uniform initial temp of sphere, K
Tinf = 773      # surrounding fluid or gas temp, K
tmax = 1.0      # max time, s
h = 375         # heat transfer coefficent, W/(m2 K)

# Initial Calculations
#------------------------------------------------------------------------------

ro = (d/2)                          # radius of sphere (a.k.a outer radius), m
rs = ro/ro                          # dimensionless surface radius, (-)
rc = 1e-12/ro                       # dimensionless center radius, (-)

alpha = kw/(rhow*cpw)               # thermal diffusivity biomass, m^2/s
t = np.arange(0, tmax+0.002, 0.002) # time range for simulation, s
z = np.arange(0, 1250, 0.1)         # range to evaluate the zeta, Bi equation

Bi = (h*ro)/kw                      # Biot number, (-)
Fo = (alpha * t) / (ro**2)          # Fourier number, (-)

# Sphere Temperature Profiles
#------------------------------------------------------------------------------

b = 2   # shape factor where 2 sphere, 1 cylinder, 0 slab

# surface temperature where ro for outer surface
thetaRo = theta(rs, b, z, Bi, Fo)   # dimensionless temperature profile
T_o = Tinf + thetaRo*(Ti-Tinf)      # convert theta to temperature in Kelvin, K

# center temperature where r for center
thetaR = theta(rc, b, z, Bi, Fo)    # dimensionless temperature profile
T_r = Tinf + thetaR*(Ti-Tinf)       # convert theta to temperature in Kelvin, K

# Cylinder Temperature Profiles
#------------------------------------------------------------------------------

b = 1   # shape factor where 2 sphere, 1 cylinder, 0 slab

# surface temperature where ro for outer surface
thetaRo = theta(rs, b, z, Bi, Fo)   # dimensionless temperature profile
To_cyl = Tinf + thetaRo*(Ti-Tinf)   # convert theta to temperature in Kelvin, K

# center temperature where r for center
thetaR = theta(rc, b, z, Bi, Fo)    # dimensionless temperature profile
Tr_cyl = Tinf + thetaR*(Ti-Tinf)    # convert theta to temperature in Kelvin, K

# Slab Temperature Profile
#------------------------------------------------------------------------------

b = 0   # shape factor where 2 sphere, 1 cylinder, 0 slab

# surface temperature where ro for outer surface
thetaRo = theta(rs, b, z, Bi, Fo)   # dimensionless temperature profile
To_slab = Tinf + thetaRo*(Ti-Tinf)  # convert theta to temperature in Kelvin, K

# center temperature where r for center
thetaR = theta(rc, b, z, Bi, Fo)    # dimensionless temperature profile
Tr_slab = Tinf + thetaR*(Ti-Tinf)   # convert theta to temperature in Kelvin, K

# Plot Results
#------------------------------------------------------------------------------

# configure y-axis based on cooling or heating simulation
if Ti > Tinf:
    # for a cooling process where Ti=773K and Tinf=300K
    Th = Ti
    ylim =[Ti+20, Tinf-20]
else:
    # for a heating process where Ti=300K and Tinf=773K
    Th = Tinf
    ylimRange = [Ti-20, Tinf+20]
    
py.figure(1)
py.plot(t, T_o, '-r', lw=2, label='surface')
py.plot(t, T_r, '--r', lw=2, label='center')
py.title('Sphere')
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

py.figure(2)
py.plot(t, To_cyl, '-b', lw=2, label='surface')
py.plot(t, Tr_cyl, '--b', lw=2, label='center')
py.title('Cylinder')
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

py.figure(3)
py.plot(t, To_slab, '-g', lw=2, label='surface')
py.plot(t, Tr_slab, '--g', lw=2, label='center')
py.title('Slab')
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
