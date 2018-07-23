#!/usr/bin/python3

import sys
import spiceypy as spice
import numpy as np
import scipy.special as scs
import datetime as dt
import calendar
import matplotlib.pyplot as plt

# year of the tidal forcing
year = int(sys.argv[1])

# path to output directory
path = "input/"

# print toolkit version numbers
print(spice.tkvrsn('TOOLKIT'))

# load kernels
spice.furnsh("meta_kernel")

def lon_lat_r(body, time, earth_radius=6371e3):

    """
    Use SPICE software to get longitude and latitude of the point on Earth at
    which <body> is in zenith and the distance between the barycenters of Earth
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
    Use SPICE software to get the product of the gravitational constant and the
    mass of <body>. Input is SPICE name, e.g. "SUN" or "MOON".
    """
    _, GM = spice.bodvrd(body, "GM", 1)
    return GM[0]*1e9 # convert km to m

def mu(lon1, lat1, lon2, lat2):
    """
    Calculate the cosine of the zenith angle alpha (using the spherical law of
    cosines).
    """
    return np.sin(lat1)*np.sin(lat2)+np.cos(lat1)*np.cos(lat2)*np.cos(lon2-lon1)

def tidal_potential(lon, lat, time, earth_radius=6371e3, h2=0.61, k2=0.30):
    """
    Calculate the tidal potential from the orbital position of the sun and moon
    (using the second-order Legendre polynomial only, see e.g. Munk and
    Cartwright, 1966, Philos. T. Roy. Soc. A), corrected for the solid earth
    tide using the Love numbers h2 and k2 (e.g. Cartwright, 1977, Rep. Prog.
    Phys.).

    Input:
      lon, lat     - longitude and latitude at which to compute the tidal
                     potential (in radians)
      time         - date and time in datetime format (in UTC)
    Optional input:
      earth_radius - Earth's radius (default 6371e3 m)
      h2, k2       - Love numbers (default h2 = 0.61, k2 = 0.30)
    Output:
      V_tide       - tidal potential
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
    sun_V = -sun_GM/sun_r*sun_xi**2*scs.legendre(2)(sun_mu)
    moon_V = -moon_GM/moon_r*moon_xi**2*scs.legendre(2)(moon_mu)

    # make solid earth tide correction
    sun_V = (1+k2-h2)*sun_V
    moon_V = (1+k2-h2)*moon_V

    return sun_V + moon_V

# number of days in the specified year
days = 366 if calendar.isleap(year) else 365

# map to output file (overwrites existing file!)
V_out = np.memmap("{:s}tide_{:4d}".format(path, year), dtype=np.dtype(">f"),
        mode="w+", shape=(days*24,181,360))

# loop over hours
for i in range(days*24):

    # time
    time = dt.datetime(year,1,1,0,0,0) + dt.timedelta(hours=i)
    print(time)

    # lon/lat
    lon = np.deg2rad(np.arange(0, 360))[None,:]
    lat = np.deg2rad(np.arange(-90, 91))[:,None]

    # get tidal forcing and save to file (using memmap)
    V_out[i,:,:] = tidal_potential(lon, lat, time)

    # save figure
    plt.imsave("fig/pot_{:4d}_{:05d}.png".format(year, i), V_out[i,:,:],
            vmin=-5, vmax=5, dpi=300)
