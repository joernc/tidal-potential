# This script adds tidal forcing in the form of surface pressure to the JRA
# fields of atmospheric surface pressure. It calculates the tidal potential from
# knowledge of the sun and moon at a given moment, and it corrects for the
# solid-earth tide. See, for example, Wunsch: Modern Observational Physical
# Oceanography, 2015, chapter 6 for the relevant theory.

# This script makes use of the NASA software SPICE. It thus requires the package
# spiceypy (github.com/AndrewAnnex/SpiceyPy), which is a Python wrapper to that
# software. Also required are the kernels listed in <metakernel.txt>, which can
# be downloaded from naif.jpl.nasa.gov/pub/naif/generic_kernels/.

import spiceypy as spice
import numpy as np
import scipy.special as scs
import datetime as dt
import calendar
import matplotlib.pyplot as plt

# path to JRA-55 files (also used as output folder!)
path = "/nobackup1/joernc/patches/jra55/"

# year of the JRA-55 file tidal forcing should be added to
year = 2014

# reference density of the ocean simulation
rho = 1029.

# print toolkit version numbers
print spice.tkvrsn('TOOLKIT')

# load kernels
spice.furnsh("meta_kernel.txt")

def lon_lat_r(body, time, earth_radius=6371e3):

    """
    Use SPICE software to get longitude and latitude of the point on Earth at
    which <body> is in zenith, and the distance between the barycenters of Earth
    and <body> at a specified <time> (in UTC). We're pretending the earth is a
    sphere here, as is commonly done in tidal studies (and presumably is done in
    the ocean model this is fed into).

    Input:
      body         - SPICE name of the body (e.g. "SUN" or "MOON")
      time         - date and time in datetime format (UTC)
      earth_radius - radius of the earth (default 6371e3 m)
    Output:
      lon, lat     - longitude and latitude at which <body> is in zenith
      r            - range of the body (distance between barycenters)
    """

    # convert UTC time to ephemeris time
    time = spice.str2et(str(time) + " UTC")

    # get position in rectangular coordinate system (body frame of Earth)
    pos_rec, _ = spice.spkpos(body, time, "ITRF93", "NONE", "EARTH")

    # transform to geodetic coordinates assuming zero flatness (to get lat/lon)
    pos_geo = spice.recgeo(pos_rec, earth_radius/1e3, 0)

    # transform to range, right ascension, declination (to get range)
    pos_rad = spice.recrad(pos_rec)

    # isolate coordinates we'll need (and convert from km to m)
    lon = pos_geo[0]
    lat = pos_geo[1]
    r = pos_rad[0]*1e3

    return lon, lat, r

def GM(body):
    """
    Use SPICE software to get the product of the gravitational constant and thei
    mass of <body>. Input is SPICE name, e.g. "SUN" or "MOON".
    """
    _, GM = spice.bodvrd(body, "GM", 1)
    return GM[0]*1e9 # convert km to m

def mu(lon1, lat1, lon2, lat2):
    """
    Calculate the cosine of the zenith angle alpha (using the spherical law of
    cosines.
    """
    return np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1)

def pres_tide(lon, lat, time, rho, earth_radius=6371e3, h2=0.61, k2=0.30):
    """
    Calculate the surface pressure representing the body force due to the tidal
    potential. The tidal potential is calculated from the orbital position of
    the sun and moon (using the second-order Legendre polynomial only, see e.g.
    Munk and Cartwright, 1966, Philos. T. Roy. Soc. A), corrected for the solid
    earth tide using the Love numbers h2 and k2 (e.g. Cartwright, 1977, Rep.
    Prog. Phys.), and converted into a pressure by multiplication by rho (the
    reference density in the ocean simulation that is to be forced with this).

    Input:
      lon, lat     - longitude and latitude at which to compute pressure (in
                     radians)
      time         - date and time in datetime format (in UTC)
      rho          - reference density of ocean simulation
    Optional input:
      earth_radius - Earth's radius (default 6371e3 m)
      h2, k2       - Love numbers (default h2 = 0.61, k2 = 0.30)
    Output:
      p_tide       - surface pressure field representing tidal forcing
    """

    # get position in rectangular coordinate system
    sun_lon, sun_lat, sun_r = lon_lat_r("SUN", time)
    moon_lon, moon_lat, moon_r = lon_lat_r("MOON", time)

    # cosine of zenith angle alpha
    sun_mu = mu(lon, lat, sun_lon, sun_lat)
    moon_mu = mu(lon, lat, moon_lon, moon_lat)

    # parallaxes
    sun_xi = earth_radius/sun_r
    moon_xi = earth_radius/moon_r

    # GM
    sun_GM = GM("SUN")
    moon_GM = GM("MOON")

    # potentials
    sun_V = sun_GM/sun_r*sun_xi**2*scs.legendre(2)(sun_mu)
    moon_V = moon_GM/moon_r*moon_xi**2*scs.legendre(2)(moon_mu)

    # convert to pressure and make solid earth tide correction
    sun_p = (1+k2-h2)*rho*sun_V
    moon_p = (1+k2-h2)*rho*moon_V

    return sun_p + moon_p

def pres_jra(time, path):
    """
    Load JRA-55 surface pressure field located in <path> and interpolate to the
    time given by <time>. This is assuming <time> is given as a datetime object
    and is of the form hh:30:00.
    """

    # decide whether and how to interpolate
    if time.hour % 3 == 0: # one hour before JRA time step
        # times and weights of JRA to be interpolated between
        t0 = time - dt.timedelta(hours=2)
        t1 = time + dt.timedelta(hours=1)
        w0 = 1./3
        w1 = 2./3
    elif time.hour % 3 == 1: # matches JRA time stamp
        # times and weights of JRA (no interpolation needed here)
        t0 = time
        t1 = time
        w0 = 1
        w1 = 0
    elif time.hour % 3 == 2: # one hour after JRA time step
        # times and weights of JRA to be interpolated between
        t0 = time - dt.timedelta(hours=1)
        t1 = time + dt.timedelta(hours=2)
        w0 = 2./3
        w1 = 1./3

    # map to atmospheric surface pressure binary files
    days0 = 366 if calendar.isleap(t0.year) else 365
    days1 = 366 if calendar.isleap(t1.year) else 365
    pres0 = np.memmap("{:s}jra55_pres_{:4d}".format(path, t0.year),
            dtype=np.dtype(">f"), mode="r", shape=(days0*24/3,320,640))
    pres1 = np.memmap("{:s}jra55_pres_{:4d}".format(path, t1.year),
            dtype=np.dtype(">f"), mode="r", shape=(days1*24/3,320,640))

    # interpolate pressure field
    delt0 = t0 - dt.datetime(t0.year,1,1,0,30,0)
    delt1 = t1 - dt.datetime(t1.year,1,1,0,30,0)
    i0 = delt0.days*24 + delt0.seconds/3600
    i1 = delt1.days*24 + delt1.seconds/3600
    pres = w0*pres0[i0/3,:,:] + w1*pres1[i1/3,:,:]

    return pres

# load JRA-55 grid
f = np.load("jra55_grid.npz")
lat = np.deg2rad(f["lat"])
lon = np.deg2rad(f["lon"])

# number of days in the specified year
days = 366 if calendar.isleap(year) else 365

# map to output file (overwrites existing file, change to "r+"?)
p_out = np.memmap("{:s}jra55_pres_tide_{:4d}".format(path, year),
        dtype=np.dtype(">f"), mode="w+", shape=(days*24,320,640))

# loop over hours
for i in range(days*24):

    # time
    time = dt.datetime(year,1,1,0,30,0) + dt.timedelta(hours=i)
    print time

    # get tidal forcing
    p_tide = pres_tide(lon, lat, time, rho)

    # get JRA pressure
    p_atms = pres_jra(time, path)

    # save to file (using memmap)
    p_out[i,:,:] = p_atms + p_tide

    # save figure (origin bug!)
    plt.imsave("fig/pot_{:4d}_{:05d}.png".format(year, i),
            (p_atms+p_tide)[::-1,:], vmin=900e2, vmax=1080e2, dpi=300)
