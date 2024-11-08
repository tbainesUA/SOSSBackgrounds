import os

from pathlib import Path
from datetime import datetime, timezone


import numpy as np
import pandas as pd

from astropy.io import fits
from astorpy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.time import Time

from jwst import datamodels

from skimage.restoration import denoise_nl_means, estimate_sigma
from maskfill import maskfill


from preprocessing.flatfield import apply_flat_field
from source_detection import detect_source_mask

nonref = (slice(-4, 4), slice(4, -4))


def get_current_utc():
    """Get the current UTC time formatted as a string."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]


def convert_to_phase(times: Time):
    """phase folded dates over 1 year"""
    return ((times.jd - 2451545.0) / 365.25) % 1


def convert_to_years(times: Time):
    """convert"""
    return 2000 + (times.jd - 2451545.0) / 365.25


def convert_to_day_of_year(times: Time):
    decimalyear = times.decimalyear
    return (365.25 * (decimalyear - 2022.0)) % 365.25


def apply_non_local_means(image, patch_size=5, patch_distance=6, h=0.8):
    sigma_est = np.mean(estimate_sigma(image))
    denoised_image = denoise_nl_means(
        image,
        h=h * sigma_est,
        sigma=sigma_est,
        patch_size=patch_size,
        patch_distance=patch_distance,
        preserve_range=True,
    )
    return denoised_image


def apply_denoise_nl_means(image):
    sigma_est = estimate_sigma(image)
    return denoise_nl_means(
        image,
        patch_size=5,
        patch_distance=6,
        fast_mode=True,
        h=0.8 * sigma_est,
        sigma=sigma_est,
        preserve_range=True,
    )


def apply_maskfill(image, mask=None):
    if mask is None:
        mask = np.isnan(image)
    return maskfill(image, mask=mask, smooth=False)[0]


def group_sequential_observations(df, gap_days=1):
    """
    Groups sequential observations based on a specified gap in days.

    Parameters:
        df (pd.DataFrame): DataFrame containing the observations.
        date_col (str): The column name with the observation dates in ISO format.
        gap_days (int): Number of days to define a gap. Observations with gaps greater than this will start a new group.

    Returns:
        pd.DataFrame: DataFrame with an additional 'group' column indicating sequential groups.
    """
    # Convert the date column to datetime if not already
    times = pd.to_datetime(df)

    # Calculate the time difference between each row and the previous row
    time_eff = times.diff()

    # Define the gap threshold
    gap_threshold = pd.Timedelta(days=gap_days)

    # Create a group identifier for each sequential set of observations
    group = (time_eff > gap_threshold).cumsum()

    # Drop 'time_diff' if no longer needed
    # df = df.drop(columns="time_diff")
    return group


# environment variables
zodical_csv = os.getenv("ZODICAL_CSV")
flatfield_dir = Path(os.getenv("FLATFIELD_DIR"))
sossbkg_dir = Path(os.getenv("SOSSBKG_DIR"))


# load data
zodical_stage1_df = pd.read_csv(zodical_csv)

# Ensure the output directory exists
flatfield_dir.mkdir(parents=True, exist_ok=True)

# Update the DataFrame with the processed file paths
fits_files = zodical_stage1_df["FILEPATH"]

# apply flatield step
flatfield_fits_files = apply_flat_field(fits_files, flatfield_dir, parallel=True)

# add files to the dataframe
zodical_stage1_df["FLATFIELD_FILEPATHS"] = pd.Series(flatfield_fits_files)


# Step for running through files to group them to make bkg images

# keys which to group data by
group_keys = ["PROGRAM", "TARGPROP", "OBSERVTN"]

include_pwcpos = False

if include_pwcpos:
    group_keys += ["PWCPOS"]


sossbkg_dir.mkdir(exist_ok=True, parents=True)

# grouped dataframes
grouped_dfs = zodical_stage1_df.groupby(group_keys)

# Main loop to iterate over the grouped data frames by the Program ID, target,
# and observation number. Datasets are then given a group ID based on
# observations that were taking sequentially. There are a total of 10 exposures
# per observation for the Calibration acitivty programs that should have been
# taken one after another. However, im guessing due to scheduling windows
# some data was taken a few or many days apart. We know that the background
# intensity changes over the course of the year. In addition, it's hypothesized
# that the background may vary spatailly given the PWCPOS keyword.
for keys, df in grouped_dfs:
    df["obs_groups"] = df["DATE-BEG"].transform(
        lambda x: group_sequential_observations(x)
    )

    # iterate over the observation groups to process and combine the file datasets
    for _, files in df.groupby(["obs_groups"])["FLATFIELD_FILEPATH"]:

        # load the data
        input_models = []
        for file in files:
            with fits.open(file) as hdul:
                input_models.append(datamodels.CubeModel(hdul))

        # make initial background image
        bkg = sigma_clipped_stats([dm.data[0] for dm in input_models], sigma=1, axis=0)

        # source masking
        shape = (len(files), 2040, 2040)
        source_masks = np.zeros(shape, dtype=bool)

        for i, dm in enumerate(input_models):
            data = dm.data[0][nonref]
            err = dm.err[0][nonref]
            mask = dm.dq[0][nonref] > 0
            source_masks[i] = detect_source_mask(
                data,
                err,
                mask,
                bkg=bkg[nonref],
            )

        # apply source mask and reconstruct the median clipped combined image
        bkg = sigma_clipped_stats(
            [dm.data[0][nonref] for dm in input_models],
            mask=source_masks,
            sigma=1,
            axis=0,
        )

        # apply maskfill to fill nans
        bkg = apply_maskfill(bkg, mask=np.isnan(mask))

        # apply the non local mean denoising algorithm to smooth image
        bkg_denoise = apply_denoise_nl_means(bkg)

        # FITS formating
        bkg_fits = fits.ImageHDU(data=bkg, name="BKG")
        bkg_fits.header["DESCRIP"] = "Flatfielded Empirical SOSS Background"
        bkg_dn_fits = fits.ImageHDU(data=bkg_denoise, name="BKGDN")
        bkg_dn_fits.header["DESCRIP"] = (
            "Flatfielded Empirical SOSS Background NLM denoised"
        )

        # RA and Dec pointing averages associated with the background obs
        skycoords = np.array(
            [
                (
                    fits.getval(file, "RA_V1", extname="SCI"),
                    fits.getval(file, "DEC_V1", extname="SCI"),
                )
                for file in files
            ]
        )
        skycoords_mean = np.mean(skycoords, axis=0)
        ra_obs, dec_obs = skycoords_mean
        sky_coords = SkyCoord(ra_obs, dec_obs, format="icrs", unit="deg")

        # convert to ecliptic coordinates
        ecliptic_coords = sky_coords.transform_to("geocentricmeanecliptic")
        lon = ecliptic_coords.lon.value
        lat = ecliptic_coords.lat.value

        # # average/middle time of the background

        mid_time_mjd = np.array(
            [model.meta.exposure.mid_time_mjd for model in input_models]
        )
        mid_time_mjd = Time(mid_time_mjd, format="mjd")
        day_of_year = convert_to_day_of_year(mid_time_mjd)

        outfile = f"nis-clear-gr700xd-ECLIPTIC-LON{lon:.2f}-LAT{lat:.2f}_{keys[-1]}_skybkg.fits"

        primary = fits.PrimaryHDU()

        date_created = get_current_utc()

        primary.header["DATE"] = date_created, "Date file created"
        primary.header["FILENAME"] = outfile
        primary.header["AUTHOR"] = "Tyler Baines"
        primary.header["DESCRIP"] = "Empirical SOSS Background"
        primary.header["FLATCORR"] = True, "Flatfield corrected images"
        primary.header["NSAMPS"] = (
            len(files),
            "Number of images for empirical background",
        )
        primary.header["DAY"] = np.median(day_of_year), "median day of year"
        primary.header["MINDAY"] = np.min(day_of_year), "min day of year"
        primary.header["MAXDAY"] = np.max(day_of_year), "max day of year"

        outfits = fits.HDUList([primary, bkg_fits, bkg_denoise])

        outfits.writeto(sossbkg_dir / outfile, overwrite=True)
