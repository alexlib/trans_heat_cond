# Implicit numerical solution of 1D transient heat conduction
# b = 1 cylinder, b = 2 sphere
# convection at surface & no heat of reaction
# see Ozisik1993, Ch.12, pg.459
# includes surface plot & radius plot

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# modules
import numpy as np
import matplotlib.pyplot as py
from mpl_toolkits.mplot3d.axes3d import Axes3D

import time

ti = time.clock()   # start time

# Parameters from Papadikis2010a Table 1
# -------------------------------------------------------------------------
rho = 700       # density of wood, kg/m^3
d = 0.035e-2    # wood particle diameter, m
cpw = 1500      # biomass specific heat capacity, J/kg*K
kw = 0.105      # biomass thermal conductivity, W/m*K
h = 375         # heat transfer coefficient, W/m^2*K
Ti = 300        # initial particle temp, K
Tinf = 773      # ambient temp, K

# Implicit Numerical Model where b = 1 cylinder and b = 2 sphere
# -------------------------------------------------------------------------
nt = 499       # number of time steps
tmax = 0.8      # max time, s
dt = tmax/nt    # time step, s
t = np.arange(0,tmax+dt,dt)

nr = 49     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

b = 2   # run model as a cylinder (b = 1) or as a sphere (b = 2)

alpha = kw/(rho*cpw)    # thermal diffusivity, alfa = kw / rho*cp, m^2/s
Fo = alpha*dt/(dr**2)   # Fourier number, Fo = alfa*dt / dr^2, (-)
Bi = h*dr/kw            # Biot numbmer, Bi = h*dr / kw, (-)

# create array [TT] to store temperature values
# note that row = time step, column = node
TT = np.zeros((len(t), m))
TT[0, :] = Ti   # first row is initial temperature of sphere or cylinder

# build tridiagonal coefficient matrix [A] and initial column vector {C}
A = np.zeros((m,m)) # pre-allocate [A] array
C = np.zeros((m,1)) # pre-allocate {C} vector

A[0, 0] = 1 + 2*(1+b)*Fo
A[0, 1] = -2*(1+b)*Fo
C[0, 0] = Ti

for i in range(1, m-1):
    A[i, i-1] = -Fo*(1 - b/(2*(i+1)))   # Tm-1
    A[i, i] = 1 + 2*Fo                  # Tm
    A[i, i+1] = -Fo*(1 + b/(2*(i+1)))   # Tm+1
    C[i, 0] = Ti

A[m-1, m-2] = -2*Fo
A[m-1, m-1] = 1 + 2*Fo*(1 + Bi + (b/(2*m))*Bi)
C[m-1, 0] = Ti + 2*Fo*Bi*(1 + b/(2*m))*Tinf

# solve system of equations [A]{T} = {C} for column vector {T}
for i in range(1, nt+1):
    T = np.linalg.solve(A,C)
    C = T.copy()
    C[m-1, 0] = T[m-1, 0] + 2*Fo*Bi*(1 + b/(2*m))*Tinf
    TT[i, :] = T.T

# Plot results
# -------------------------------------------------------------------------
rr = np.linspace(0, r, m)   # radius vector
rn = rr/r                   # normalized radius vector

# surface & center temp. vs time plot
py.figure(1)
py.plot(t, TT[:, m-1], '-k', label='surface')
py.plot(t, TT[:, 0], '--k', label='center')
py.axhline(Tinf, color='r', linestyle='-.')
py.ylabel('Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.ylim([Ti-20, Tinf+20])
py.xlim([0, tmax])
py.grid()
py.show()

# radius vs temp. at different times
t1 = len(t)*(1/5)
t2 = len(t)*(2/5)
t3 = len(t)*(3/5)
t4 = len(t)*(4/5)
t5 = len(t)*(5/5)
print('times ', t1*dt, t2*dt, t3*dt, t4*dt, t5*dt)

py.figure(2)
py.plot(rn, TT[t5-1, :], label='%.2f s' % t[t5-1])
py.plot(rn, TT[t4, :], label='%.2f s' % t[t4])
py.plot(rn, TT[t3, :], label='%.2f s' % t[t3])
py.plot(rn, TT[t2, :], label='%.2f s' % t[t2])
py.plot(rn, TT[t1, :], label='%.2f s' % t[t1])
py.ylabel('Temperature (K)')
py.xlabel('Radius (-)')
py.legend(loc='best', numpoints=1)
py.show()

# radius, time, temp. surface plot 
rn, tt = np.meshgrid(rn, t)

fig = py.figure(3)
ax = fig.add_subplot(1, 1, 1, projection='3d')
p = ax.plot_surface(rn, tt, TT, rstride=10, cstride=10, cmap='RdYlGn_r', linewidth=0, antialiased=False)
ax.set_xlabel('Radius (-)')
ax.set_ylabel('Time (s)')
ax.set_zlabel('Temp (K)')
py.show()

# elapsed time for entire file
print('time', time.clock()-ti, 'seconds')