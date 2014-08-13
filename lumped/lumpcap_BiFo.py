# Lumped Capacitance Model for Sphere using Biot and Fourier numbers
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

Ti = 300    # initial temperature, K
Tinf = 773  # ambient temperature, K

t = np.linspace(0, 3)   # time range, s

# Calculations for Sphere
# -----------------------------------------------------------------------------

A = np.pi*(d**2)        # surface area sphere pi*D^2, m^2
V = np.pi*(d**3)/6      # volume of sphere pi*D^3/6, m^3

Lc = V/A          # characteristic length for sphere, m
Bi = (h*Lc)/k     # Biot number Eq 5.10

alpha = k/(rho*c)           # thermal diffusivity, m^2/s
Fo = (alpha*t)/(Lc**2) # Fourier number Eq 5.12

phi = np.exp(-Bi*Fo)    # dimensionless temperature, Eq 5.6
T = Tinf+(Ti-Tinf)*phi  # temperature, K

# Plot Results
# -----------------------------------------------------------------------------

py.figure(3)
py.plot(t, phi)
py.title('Sphere')
py.ylabel(r'$\Theta^*$', rotation=0)
py.xlabel('t (s)')
py.grid()
py.show()

py.figure(4)
py.plot(t, T)
py.axhline(y=773, color='k', linestyle='--')
py.ylim(200,800)
py.title('Sphere')
py.ylabel('T (K)', rotation=0)
py.xlabel('t (s)')
py.grid()
py.show()