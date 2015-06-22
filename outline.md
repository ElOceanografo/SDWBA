## SDWBA package ##

The goal of this package is to provide a set of easy-to-use tools to calculate the acoustic backscattering cross-sections of marine organisms using the (stochastic) distorted-wave Born approximation (the DWBA or SWDBA).  These models discretize zooplankton or other scatterers as a series of cylindrical sections, and are efficient and accurate for fluid-like organisms, including krill and copepods.

Methods, parameters, and modifications for these calculations have been previously published in a number of different papers.  However, the code necessary to perform these calculations in practice is *not* always published.  When it is, may be written for specific organisms or geometries, making it difficult to generalize and extend.  Updates to methods published in the literature may not find their way into practice as easily as they would if they could be added to a standard set of 

We aim to provide functions which calculate the backscattering cross-section (or target strength) of a fluid-like scatterer with arbitrary shape, material properties, and orientation.  These will be made available in several popular languages for technical computing, including at least R, Matlab, and Python, with tests to ensure identical results in each.  We will also include functions for resizing animals and resampling their shapes to use an adequate number of scattering elements. The code will be clearly implemented, organized, and well-documented, and referenced to the original papers.


### Structures/classes ###

Though the details will vary depending on the language, each implementation will include a data structure specifying the body shape and material properties of the scattering organism.  Basically, this will contain the following vectors which digitize the critter at a series of points along its body:

* x, y, z : points along the animal's centerline.  These are defined in the animal's "standard"
      coordinates, e.g. with the animal's main axis horizontal and the tail at (0, 0, 0).
* a : radii corresponding to those points
* density : density at each point in the animal's body
* soundspeed : sound speed at each point in the animal's body

### Functions ###

The following functions will be available in each language (though Python syntax is shown).  Arguments are briefly described below each function.

---

`load_scatterer(filename, colnames)`

Reads in the shape and properties of a scatterer from a .csv file on disk.  Returns a `Scatterer` object/data structure.

* `filename` : path to .csv file with scatterer shape and material properties
* `colnames` : (optional) dictionary saying which column is which, in case their names are different from those expected in the other code.

---

`save_scatterer(scatterer, filename)`

Saves the shape and properties of a scatterer to a .csv file on disk.

* `scatterer` : a `Scatterer` object.
* `filename` : path to .csv file with scatterer shape and material properties

---

`resample_scatterer(scatterer, n=None, s=None, freq=None)`

Resamples the digitization points along the scatterer's centerline.

* `scatterer` : a `Scatterer` object.
* `n` : The desired number of points.  If not supplied, it will be calculated as the minimum number necessary to represent the animal based on its length and the specified frequency.
* `s` : A vector of digitization locations, as lengths along the animal's centerline.  Can be supplied instead of `n`.
* `freq` : Acoustic frequency in Hz.  Only used if `n` or `s` is not provided.

---

`resize.scatterer(scatterer, strech, resample=True, freq=None)`

Resize a scatterer by a constant factor.  Returns a `Scatterer` object with the new shape.

* `scatterer` : a `Scatterer` object.
* `stretch` : Either a scalar specifying a constant factor, or a length-2 vector specifying how much to stretch the scatterer in length and width.
* `resample` : Should the digitizing points be recalculated if the new length makes them too sparse with respect to the acoustic wavelength?
* `freq` : Acoustic frequency in Hz.  Only used if `resample=True`.

----

`rotate(scatterer, tilt=0, roll=0, yaw=0)`

Rotate a scatterer in space, starting from its initial coordinate system.  Returns a new `Scatterer` object.

* `scatterer` : a `Scatterer` object.
* `tilt`, `roll`, `yaw` : angles (in degrees) for rotating the scatterer.  These are applied in order.

---

`backscatter.xsection(scatterer, k, sound_speed, density, phase_sd=0, TS=False)`

Calculate the backscattering cross-section (i.e. $\sigma_{bs}$), in $m^2$, of a scatterer with the incident acoustic field described by `k`.

* `scatterer` : a `Scatterer` object
* `k` : the length-3 vector of the acoustic wavenumber. For a downward-looking echosounder in a normal z-up coordinate system, this is `[0, 0, -2 * pi * f]` (where `f` is the frequency) 
* `soundspeed` : sound speed in water, in m/s.
* `density` : density of surrounding water, in kg/m^3
* `phase_sd` : standard deviation of the phase of the random scattering component, in radians.  Defaults to zero.
* `TS` : Should the result be returned in logarithmic form, i.e. as dB referenced to 1 m^2? Defaults to False.

---

`frequency.scattering(scatterer, freq_start, freq_stop, step=1e3, *args, **kwargs)`

Calculate the backscattering cross-section of the scatterer across a range of frequencies.
* `scatterer` : a `Scatterer` object
* `freq_start`, `freq_stop` : minimum and maximum frequencies desired
* `step` : Frequency interval.  Defaults to 1 kHz.
* `*args, **kwargs` : Additional arguments passed to `backscatter.xsection`.

---

`tilt.scattering(scatterer, tilt_start, tilt_stop, step=1, *args, **kwargs)`

Calculate $\sigma_{bs}$ for a scatterer across a range of tilt angles (i.e., departures from "horizontal", however that is naturally defined for the critter).

* `scatterer` : a `Scatterer` object
* `tilt_start`, `tilt_stop` : minimum and maximum tilt angles (in degrees) desired
* `step` : Angle interval.  Defaults to 1 degree.* 
* `*args, **kwargs` : Additional arguments passed to `backscatter.xsection`.

### Data ###

This package will separate the calculation of an organism's backscatter or target strength from the definition of its shape and material properties.  The advantage of this approach is that it is easy to define models for new organisms. 

The package should come "batteries included" with data files containing the shapes of some common  marine zooplankton:

* Antarctic krill
* Other krill species
* Copepods
