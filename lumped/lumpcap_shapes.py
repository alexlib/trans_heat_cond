# Lumped Capacitance Model comparing sphere, cylinder, cube shapes
# uses sphere diameter to determine cylinder and cube of same volume as sphere
# Incropera2011, Ch.5 Transient Conduction, pg.280-286

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Parameters
# -----------------------------------------------------------------------------

h = 375     # convection heat transfer coefficient, W/m^2*K
c = 1500    # specific heat, J/kg*K
k = 0.20    # thermal conductivity, W/m*K
rho = 700   # density, kg/m^3

d = 350e-6  # particle diameter, m (e-6 for microns)
H = 2/3*d                   # cylinder height of equal volume
a = ((np.pi*d**3)/6)**(1/3) # cube side of equal volume

Ti = 300    # initial temperature, K
Tinf = 773  # ambient temperature, K

t = np.linspace(0, 2)   # time range, s

# Functions
# -----------------------------------------------------------------------------

# returns the volume of a sphere, cylinder, and cube
def volume(d, h, a):
    sphere = (np.pi*d**3)/6
    cylinder = (np.pi*d**2)/4*h
    cube = a**3
    return sphere, cylinder, cube
    
    
# returns the surface area of a sphere, cylinder, and cube
def surfArea(d, h, a):
    sphere = np.pi*d**2
    cylinder = np.pi*d*h + (np.pi*d**2)/2
    cube = 6*a**2
    return sphere, cylinder, cube
    
    
# return the characteristic length of a sphere, cylinder, cube
def charLength(Vsph, Vcyl, Vcube, Asph, Acyl, Acube):
    sphere = Vsph/Asph
    cylinder = Vcyl/Acyl
    cube = Vcube/Acube
    return sphere, cylinder, cube
    
    
# Calculations
# -----------------------------------------------------------------------------

# volume, surface area, and characteristic length
Vsph, Vcyl, Vcube = volume(d, H, a)
Asph, Acyl, Acube = surfArea(d, H, a)
Lsph, Lcyl, Lcube = charLength(Vsph, Vcyl, Vcube, Asph, Acyl, Acube)

# Biot number
Bi_sph = (h*Lsph)/k
Bi_cyl = (h*Lcyl)/k
Bi_cube = (h*Lcube)/k

# thermal diffusivity, m^2/s
alpha = k/(rho*c)

# Fourier number
Fo_sph = (alpha*t)/(Lsph**2)
Fo_cyl = (alpha*t)/(Lcyl**2)
Fo_cube = (alpha*t)/(Lcube**2)

# dimensionless temperature
phi_sph = np.exp(-Bi_sph*Fo_sph)
phi_cyl = np.exp(-Bi_cyl*Fo_cyl)
phi_cube = np.exp(-Bi_cube*Fo_cube)

# temperature, K
T_sph = Tinf+(Ti-Tinf)*phi_sph
T_cyl = Tinf+(Ti-Tinf)*phi_cyl
T_cube = Tinf+(Ti-Tinf)*phi_cube

# Plot Results
# -----------------------------------------------------------------------------

py.figure(5)
py.plot(t, phi_sph, label='sphere')
py.plot(t, phi_cyl, label='cylinder')
py.plot(t, phi_cube, label='cube')
py.title('Dimensionless Temperature')
py.ylabel(r'$\Theta^*$', rotation=0)
py.xlabel('t (s)')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()

py.figure(6)
py.plot(t, T_sph, label='sphere')
py.plot(t, T_cyl, label='cylinder')
py.plot(t, T_cube, label='cube')
py.axhline(y=773, color='k', linestyle='--')
py.title('Temperature')
py.ylabel('T (K)', rotation=0)
py.xlabel('t (s)')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()

