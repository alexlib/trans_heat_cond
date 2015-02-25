"""
Lumped capacitance method for solid sphere using Biot and Fourier numbers.

Reference:
Bergman, Lavine, Incropera, Dewitt 2011, Ch. 5, pg. 280-286
"""
 
# Modules
# -----------------------------------------------------------------------------

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

tau = (rho*V*c)/(h*A)   # thermal time constant Eq 5.7, s

Lc = V/A                # characteristic length for sphere, m
Bi = (h*Lc)/k           # Biot number Eq 5.10, (-)

alpha = k/(rho*c)       # thermal diffusivity, m^2/s
Fo = (alpha*t)/(Lc**2)  # Fourier number Eq 5.12, (-)

phi = np.exp(-Bi*Fo)    # dimensionless temperature Eq. 5.13, (-)
T = Tinf+(Ti-Tinf)*phi  # temperature, K

# time for solid to reach some temperature Eq 5.5, s
ts = tau*np.log((Ti-Tinf)/(772-Tinf))

# Print Results
# -----------------------------------------------------------------------------

print('tau (s) = ', tau)
print('Bi (-) = ', Bi)
print('ts (s) = ', ts)

# Plot Results
# -----------------------------------------------------------------------------

py.figure(1)
py.plot(t, phi, lw=2)
py.title('Sphere')
py.ylabel(r'$\Theta$ / $\Theta_i$ (-)')
py.xlabel('t (s)')
py.grid()
py.show()

py.figure(2)
py.plot(t, T, lw=2)
py.axhline(y=773, color='k', linestyle='--')
py.ylim(200,800)
py.title('Sphere')
py.ylabel('T (K)')
py.xlabel('t (s)')
py.grid()
py.show()
