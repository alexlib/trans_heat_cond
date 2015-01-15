## Analytical Approach for 1D Transient Heat Conduction

The [Analytical Model](http://nbviewer.ipython.org/github/pyrolysis/trans_heat_cond/blob/master/analytical/analytical.ipynb) uses the analytical solution for 1D transient heat conduction in a solid sphere, cylinder, or slab shape with convection at the surface and no heat generation. Run the `analytical.py` file for plots of the temperature profiles at the center and surface of the wood particle.  

Model - `analytical.py`  
Functions - `funcTheta.py`, `funcRoots.py`, `funcZeta.py`  
Examples - `bessel.py`, `roots.py`  

### analytical.py
Plots the analytical solutions of 1D transient heat conduction in a solid sphere, cylinder, and slab with convection at the surface and no heat generation within the solid. Available as a Python file or iPython notebook. The notebook can be viewed in any web browser by clicking the "Analytical Model" link above.

### bessel.py
Example of using SciPy for calculating bessel functions of the first kind. Compare to the values listed in the Bergman 2011 book.

### funcRoots.py
Function that returns the positive roots of the zeta, Bi equation of the analytical solution for 1D transient heat conduction of a sphere, cylinder, or slab shape.

### funcTheta.py
Function to return the theta (dimensionless temperature) profile at a certain dimensionless point (r) in a 1D solid sphere, cylinder, or slab shape.

### funcZeta.py
Functions for each zeta, Bi equation for sphere, cylinder, and slab cases for use in the 1D transient heat conduction analytical solution.

### roots.py
Plot the zeta, Bi equation and its positive roots. For use in the analytical transient heat conduction of a 1D sphere, cylinder, or slab.

### References
* Bergman, T.L. et al., 2011. Fundamentals of Heat and Mass Transfer 7th ed., John Wiley and Sons, Inc.
* Recktenwald, G., 2006. Transient One-Dimensional Heat Conduction in a Convectively Cooled Sphere, pp.1–13.
* Papadikis, K., Gu, S. and Bridgwater, A.V., 2010. Computational modelling of the impact of particle size to the heat transfer coefficient between biomass particles and a fluidised bed. Fuel Processing Technology, 91(1), pp.68–79.

