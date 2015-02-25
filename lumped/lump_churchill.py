"""
Lumped capacitance method for sphere applied to Churchill correlation

References:
1) Bergman, Lavine, Incropera, Dewitt 2011, Ch. 5, pg. 280-286
2) Churchill 1974
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
k = 0.70    # thermal conductivity, W/m*K
rho = 700   # density, kg/m^3
d = 350e-6  # particle diameter, m (e-6 for microns)

Ti = 300    # initial temperature, K
Tinf = 773  # ambient temperature, K

t = np.linspace(0, 2)   # time range, s
T5 = np.zeros(len(t))  # vectors to store T{t}
T10 = np.zeros(len(t))
T100 = np.zeros(len(t))
T200 = np.zeros(len(t))

# Calculations Sphere
# -----------------------------------------------------------------------------

A = np.pi*(d**2)        # surface area sphere pi*D^2, m^2
V = np.pi*(d**3)/6      # volume of sphere pi*D^3/6, m^3

# Functions
# -----------------------------------------------------------------------------

def T(h,A,rho,c,V,n,t):
    tm = (h*A)/(rho*c*V)
    To = Tinf + (Ti-Tinf)*np.exp(-tm*t)
    return To / (1+(To/Tinf)**n)**(1/n)
    
    
# Calculations
# -----------------------------------------------------------------------------

k = 0

for i in t:
    T5[k] = T(h,A,rho,c,V,5,i)
    T10[k] = T(h,A,rho,c,V,10,i)
    T100[k] = T(h,A,rho,c,V,100,i)
    T200[k] = T(h,A,rho,c,V,1000,i)
    k+=1
    
    
# Plot results
# -----------------------------------------------------------------------------

py.figure(10)
py.plot(t, T5, '-m', label='n=5')
py.plot(t, T10, '-b', label='n=10')
py.plot(t, T100, '-g', label='n=100')
py.plot(t, T200, '-r', label='n=200')
py.axhline(y=773, color='k', linestyle='--')
py.legend(loc='best', numpoints=1)
py.title('Sphere (lumped cap churchill)')
py.xlabel('t')
py.ylabel('T')
py.grid()
py.show()
