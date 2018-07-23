The script `tidal_potential.py` calculates the tidal forcing for a given year and saves it in hourly increments. The script calculates the tidal potential from knowledge of the sun and moon at a given moment, and it corrects for the solid-earth tide. See, for example, Wunsch: Modern Observational Physical Oceanography, 2015, chapter 6 for the relevant theory.

This script makes use of the NASA software SPICE. It thus requires the package [spiceypy](https://github.com/AndrewAnnex/SpiceyPy), which is a Python wrapper to that software. Also required are the kernels listed in `metakernel.txt`, which can be downloaded from [NASA's Navigation and Ancillary Information Facility](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/).

