## Structures/classes ##

* Scatterer :Representation of the body shape and material properties of a scattering organism.
  Basically a data structure to contain the following vectors which digitize the critter:
  - x, y, z : points along the animal's centerline.  These are defined in the animal's "standard"
  coordinates, e.g. with the animal's main axis horizontal and the tail at (0, 0, 0).
  - a : radii corresponding to those points
  - density : density at each point in the animal's body
  - soundspeed : sound speed at each point in the animal's body

## Functions ##

* import.scatterer(filename, colnames)
  filename : path to csv file with scatterer shape and material properties
  colnames : (optional) list/dictionary saying which column is which, in case their names are different

* plot.scatterer

* interpolate.scatterer

* rotate(scatterer, tilt=0, roll=0, yaw=0)

* backscatter.xsection(scatterer, angle, frequency, soundspeed, density, phase_sd)
  Arguments
  ---------
  scatterer : a Scatterer object
  angle : a scalar or vector describing the angle of the animal's body relative to
    the incident acoustic waves.  If a scalar, it is interpreted as a tilt angle; if a vector,
    the three components are the tilt, roll, and yaw respectively.
  frequency : the acoustic frequency in Hz.
  soundspeed : sound speed in water, in m/s.
  density : density of surrounding water, in kg/m^3
  phase_sd : standard deviation of the phase of the random scattering component, in radians

* scattering.spectrum(scatterer, start, stop, step=1)

* angle.scattering
  
## Data ##
Should come "batteries included" with shapes for generic krill, copepods, other animals...
