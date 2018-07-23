The Python 3 script `tidal_potential.py` calculates the tidal forcing for a given year and saves it in hourly increments. The script calculates the potential from knowledge of the position of the sun and moon at a given moment, and it corrects for the solid-earth tide. See, for example, Wunsch: Modern Observational Physical Oceanography, 2015, chapter 6 for the relevant theory.

This script makes use of the NASA software [SPICE](https://naif.jpl.nasa.gov/naif/). The script requires the package [spiceypy](https://github.com/AndrewAnnex/SpiceyPy), which is a Python wrapper to SPICE. The script also requires the kernels listed in `meta_kernel`, which can be downloaded from [NASA's Navigation and Ancillary Information Facility](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/).

To generate the tidal forcing for the year 2018, execute the script with:
```bash
./tidal-potential.py 2018
```
or
```bash
python3 tidal-potential.py 2018
```
The forcing fields that can be read into the ocean model are written to `input/`.
