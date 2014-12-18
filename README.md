# trans_heat_cond

Models and functions to simulate the 1D transient heat conduction of wood particles at fast pyrolysis conditions. The code is organized into 4 groups (details below) based on the approach or method used. All files are written in Python 3. For easy installation of Python 3, use the [Anaconda](http://www.continuum.io) distribution available for free from Continuum Analytics. Please read the comments in each file to fully understand the limitations and assumptions made for each example.

## analytical
Analytical solutions for 1D transient heat conduction in a solid sphere, cylinder, or slab.

## numerical
Numerical solutions for transient heat conduction in a solid sphere or cylinder.

## lumped
Lumped capacitance method for transient conduction.

## functions
Function files developed from the numerical and lumped capacitance methods. These can be imported 
into existing models. Each function is named in 'lowercase' while examples of using the function start with an 'uppercase'.

# References
The following is a list of references used the create these models. See the comments in each file for details about the particular references related to that model.

- here
- here
- here

# License
trans_heat_cond is available under the MIT license. See the LICENSE file for more info.
