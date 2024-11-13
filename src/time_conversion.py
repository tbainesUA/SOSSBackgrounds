from astropy.time import Time


def convert_to_phase(times: Time):
    """phase folded dates over 1 year"""
    return ((times.jd - 2451545.0) / 365.25) % 1


def convert_to_years(times: Time):
    """convert"""
    return 2000 + (times.jd - 2451545.0) / 365.25


def convert_to_day_of_year(times: Time):
    decimalyear = times.decimalyear
    return (365.25 * (decimalyear - 2022.0)) % 365.25
