# Lumped Capacitance Model for Sphere
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

tm = (h*A)/(rho*V*c)    # term 1
phi = np.exp(-tm*t)     # dimensionless temperature Eq 5.6

tau = (rho*V*c)/(h*A)   # thermal time constant Eq 5.7, s
print('tau = ', tau)

T = Tinf+(Ti-Tinf)*phi  # temperature at some time, K

Lc = d/6        # characteristic length, for sphere = d/6 or r/3
Bi = h*Lc/k     # Biot number Eq 5.10
print('Bi = ', Bi)

# time for solid to reach some temperature Eq 5.5, s
ts = tau*np.log((Ti-Tinf)/(772-Tinf))
print('ts = ', ts)
print('-'*20)

# plot results
py.figure(1)
py.plot(t, phi)
py.title('Sphere')
py.ylabel(r'$\Theta$ / $\Theta_i$')
py.xlabel('t (s)')
py.grid()
py.show()

py.figure(2)
py.plot(t, T)
py.axhline(y=773, color='k', linestyle='--')
py.title('Sphere')
py.ylabel('T (K)')
py.xlabel('t (s)')
py.grid()
py.show()
