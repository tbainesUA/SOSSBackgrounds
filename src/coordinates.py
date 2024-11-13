from astropy.coordinates import SkyCoord


def get_ecliptic_equatorial_coords():
    pass


def equatorial_to_ecliptic_coords(ra, dec):
    sky_coords = SkyCoord(ra=ra, dec=dec, frame="icrs", unit="deg")
    ecliptic_coords = sky_coords.transform_to("barycentricmeanecliptic")
    return ecliptic_coords.lon.degree, ecliptic_coords.lat.degree
