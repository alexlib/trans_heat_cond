# Lumped Capacitance Model for different shapes applied to Churchill correlation
# Incropera2011, Ch.5 Transient Conduction, pg.280-286
# Churchill1974

import numpy as np
import matplotlib.pyplot as py

#py.close('all')

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

# Functions
# -----------------------------------------------------------------------------

def T(Tinf, Ti, Bi_sph, Bi_cube, n, t):
    Fo_sph = (alpha*t)/(Lsph**2)
    Tsph = Tinf + (Ti-Tinf)*np.exp(-Bi_sph*Fo_sph)
    Fo_cube = (alpha*t)/(Lcube**2)
    Tcube = Tinf + (Ti-Tinf)*np.exp(-Bi_cube*Fo_cube)
    T = Tsph / (1+(Tsph/Tcube)**n)**(1/n)
    return T
    
    
# Calculations
# -----------------------------------------------------------------------------

# thermal diffusivity, m^2/s
alpha = k/(rho*c)

# sphere shape for lower bounds
Vsph = (np.pi*d**3)/6   # volume, m^3
Asph = np.pi*d**2       # surface area, m^2
Lsph = Vsph/Asph        # characteristic length, m
Bi_sph = (h*Lsph)/k     # Biot number, ~

# cube shape for upper bounds
a = ((np.pi*d**3)/6)**(1/3) # cube side having equal volume as sphere
Vcube = a**3                # volume, m^3
Acube = 6*a**2              # surface area, m^2
Lcube = Vcube/Acube         # characteristic length, m
Bi_cube = (h*Lcube)/k       # Biot number, ~

# calculate
k = 0

for i in t:
    T5[k] = T(Tinf, Ti, Bi_sph, Bi_cube, 5, i)
    T10[k] = T(Tinf, Ti, Bi_sph, Bi_cube, 10, i)
    T100[k] = T(Tinf, Ti, Bi_sph, Bi_cube, 100, i)
    T200[k] = T(Tinf, Ti, Bi_sph, Bi_cube, 200, i)
    k+=1
    
    
# Plot results
# -----------------------------------------------------------------------------

py.figure(11)
py.plot(t, T5, '-m', label='n=5')
py.plot(t, T10, '-b', label='n=10')
py.plot(t, T100, '-g', label='n=100')
py.plot(t, T200, '-r', label='n=200')
py.axhline(y=773, color='k', linestyle='--')
py.legend(loc='best', numpoints=1)
py.title('Sphere/Cube (lumped cap churchill shapes)')
py.xlabel('t')
py.ylabel('T')
py.grid()
py.show()
