# Implicit numerical solution of 1D transient heat conduction
# b = 1 cylinder, b = 2 sphere
# convection at surface & no heat of reaction
# see Ozisik1993, Ch.12, pg.459

# modules
import numpy as np
import matplotlib.pyplot as py
import transhc as hc

# Parameters from Papadikis2010a Table 1
# -------------------------------------------------------------------------
rho = 700       # density of wood, kg/m^3
d = 0.035e-2    # wood particle diameter, m
#c = 1500      # biomass specific heat capacity, J/kg*K
k = 0.105      # biomass thermal conductivity, W/m*K
h = 375         # heat transfer coefficient, W/m^2*K
Ti = 300        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Implicit Numerical Model where b = 1 cylinder and b = 2 sphere
# -------------------------------------------------------------------------
nt = 1000       # number of time steps
tmax = 0.8      # max time, s
dt = tmax/nt    # time step, s
t = np.arange(0,tmax+dt,dt)
lt = len(t)

nr = 99     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

# create array [TT] to store temperature values
# note that row = time step, column = node
T = np.zeros((lt, m))
T[0, :] = Ti   # first row is initial temperature of sphere or cylinder

# build tridiagonal coefficient matrix [A] and initial column vector {C}
A = np.zeros((m, m)) # pre-allocate [A] array
C = np.zeros((m, 1)) # pre-allocate {C} vector

# internal nodes Tm-1, Tm, Tm+1
j = np.arange(1, m-1)

# solve system of equations [A]{T} = {C} for column vector {T}
for i in range(1, nt+1):
    #c = 1112.0 + 4.85 * (T[i-1] - 273.15)
    c = np.ones(m)*1500
    A, C = hc.sphere(k, rho, c, dt, dr, h, T, Tinf, i, j, lt, m, A, C)
    TT = np.linalg.solve(A,C)
    T[i, :] = TT.T
    
    
# Plot results
# -------------------------------------------------------------------------
# surface & center temp. vs time plot
py.figure(1)
py.plot(t, T[:, m-1], '-k', label='surface')
py.plot(t, T[:, 0], '--k', label='center')
py.axhline(Tinf, color='r', linestyle='--')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.ylim([Ti-20, Tinf+20])
py.xlim([0, tmax])
py.grid()
py.show()

