## SDWBA packages ##

The goal of this package is to provide a set of easy-to-use tools to calculate the acoustic backscattering cross-sections of marine organisms using the (stochastic) distorted-wave Born approximation (the DWBA or SWDBA).  These models discretize zooplankton or other scatterers as a series of cylindrical sections, and are efficient and accurate for fluid-like organisms, including krill and copepods.

Methods, parameters, and modifications for these calculations have been previously published in a number of different papers.  However, the code necessary to perform these calculations in practice is *not* always published.  When it is, may be written for specific organisms or geometries, making it difficult to generalize and extend.  Updates to methods published in the literature may not find their way into practice as easily as they would if they could be added to a standard set of 

We aim to provide functions which calculate the backscattering cross-section (or target strength) of a fluid-like scatterer with arbitrary shape, material properties, and orientation.  These will be made available in several popular languages for technical computing, including at least R, Matlab, and Python, with tests to ensure identical results in each.  We will also include functions for resizing animals and resampling their shapes to use an adequate number of scattering elements. The code will be clearly implemented, organized, and well-documented, and referenced to the original papers.

This repository does not contain the actual code and implementations, just a description and sample datasets. The actual algorithms are available at the following locations:

* R (not fully tested): https://github.com/ElOceanografo/SDWBA.R
* Python (unfinished): https://github.com/ElOceanografo/SDWBA.py

### Data ###

This package will include files defining realistic shapes for several types of zooplankton.  At present, only Antarctic krill are supplied.
