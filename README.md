## Models and Functions for 1D Transient Heat Conduction

The folders in this repository (listed below) contain models and functions to simulate the 1D transient heat conduction of solid wood particles at fast pyrolysis conditions. The code is organized into 3 folders based on the approach or method used to calculate the heat conduction within a solid shape. Each folder contains a **README** document to provide details about the files in that particular folder.

The models and functions are written in Python 3 which is easily installed using the free [Anaconda](http://www.continuum.io) distribution provided by Continuum Analytics. This distribution includes the numerical libraries and plotting tools needed to run the models.

*Requirements: Python 3, NumPy, SciPy, and Matplotlib*

### analytical
[Analytical Model](http://nbviewer.ipython.org/github/pyrolysis/trans_heat_cond/blob/master/analytical/analytical.ipynb) - analytical solutions for 1D transient heat conduction in a solid sphere, cylinder, and slab shape. For more information, see the README file in the analytical folder or view the iPython notebook.

### numerical
Numerical Model - numerical solutions for 1D transient heat conduction in a solid sphere or cylinder.

### lumped
Lumped Model - lumped capacitance method for 1D transient heat conduction in a solid sphere, cylinder, and slab shape.

### License
Code in this repository is available under the MIT license. See the LICENSE file for more information.
