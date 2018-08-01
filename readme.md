The Python 3 script `tidal_potential.py` calculates the tidal forcing for a given year and saves it in hourly increments. The script calculates the potential from knowledge of the position of the sun and moon at a given moment, and it corrects for the solid-earth tide. See, for example, Wunsch: Modern Observational Physical Oceanography, 2015, chapter 6 for the relevant theory.

This script makes use of the NASA software [SPICE](https://naif.jpl.nasa.gov/naif/). The script requires the package [spiceypy](https://github.com/AndrewAnnex/SpiceyPy), which is a Python wrapper to SPICE. The script also requires the kernels listed in `meta_kernel`, which can be downloaded from [NASA's Navigation and Ancillary Information Facility](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/).

To generate the tidal forcing for the year 2018, execute the script with:
```bash
./tidal-potential.py 2018 IAU_EARTH
```
or
```bash
python3 tidal-potential.py 2018 IAU_EARTH
```
The forcing fields that can be read into the ocean model are written to `input/`.

SPICE uses Earth orientation data provided by a PCK kernel. The IAU_EARTH orientation model has relatively low accuracy but extends into the future indefinitely. Should higher accuracy be needed, use ITRF93 instead. To use this, be sure to add either `earth_070425_370426_predict.bpc` or `earth_latest_high_prec.bpc` to `meta_kernel`, depending on whether you need prediction or not. See the documentation of the [IAU_EARTH kernel](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc) for a discussion of its accuracy and the [documentation of the high-accuracy kernels](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/aareadme.txt) for which one to use.
