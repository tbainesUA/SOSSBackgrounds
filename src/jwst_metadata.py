import pandas as pd
from astropy.time import Time

from coordinates import equatorial_to_ecliptic_coords
from time_conversion import convert_to_day_of_year


def extract_metadata(model):
    """Model is a jwst datamodel"""

    def _extract_metadata(model, path, default):
        # Split the path and access each attribute dynamically
        value = model
        try:
            for attr in path.split("."):
                value = getattr(value, attr)

        except AttributeError:
            # Set to default if any attribute is missing
            value = default

        return value

    metadata_map = {
        "filename": ("meta.filename", None),
        "pwcpos": ("meta.instrument.pupil_position", None),
        "ra_v1": ("meta.pointing.ra_v1", None),
        "dec_v1": ("meta.pointing.dec_v1", None),
        "pa_v3": ("meta.pointing.pa_v3", None),
        "aper_pa": ("meta.aperture.position_angle", None),
        "jwst_x_bary": ("meta.ephemeris.spatial_x_bary", None),
        "jwst_y_bary": ("meta.ephemeris.spatial_y_bary", None),
        "jwst_z_bary": ("meta.ephemeris.spatial_z_bary", None),
        "jwst_x_geo": ("meta.ephemeris.spatial_x_geo", None),
        "jwst_y_geo": ("meta.ephemeris.spatial_y_geo", None),
        "jwst_z_geo": ("meta.ephemeris.spatial_z_geo", None),
        "mid_time_mjd": ("meta.exposure.mid_time_mjd", None),
    }

    metadata = {}
    for key, (path, default) in metadata_map.items():
        metadata[key] = _extract_metadata(model, path, default)
    return metadata


def parse_metadata_table(models):
    """models is a lsit of jwst data models"""

    metadata_table = pd.DataFrame(
        [extract_metadata(model) for model in models]
    )  # convert equatorial coordinates to ecliptic

    metadata_table[["ecliptic_lon", "ecliptic_lat"]] = metadata_table.apply(
        lambda x: equatorial_to_ecliptic_coords(x["ra_v1"], x["dec_v1"]),
        axis=1,
        result_type="expand",
    )

    # convert the mjd date to the day of year (0-365)
    metadata_table["day_of_year"] = metadata_table["mid_time_mjd"].apply(
        lambda x: convert_to_day_of_year(Time(x, format="mjd"))
    )

    return metadata_table
