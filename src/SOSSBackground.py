import os

from datetime import datetime, timezone
from dotenv import load_dotenv
from pathlib import Path

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from skimage.restoration import denoise_nl_means, estimate_sigma

from jwst import datamodels
from maskfill import maskfill

from flatfield import apply_flatfield
from source_detection import create_source_masks
from jwst_metadata import parse_metadata_table

load_dotenv()

nonref = (slice(4, -4), slice(4, -4))


def get_current_utc():
    """Get the current UTC time formatted as a string."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]


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


def load_fits_data(files):
    input_models = []
    for file in files:
        with fits.open(file) as hdul:
            input_models.append(datamodels.CubeModel(hdul))
    return input_models


def group_sequential_observations(df, gap_days=1):
    """
    Groups sequential observations based on a specified gap in days.

    Parameters:
        df (pd.DataFrame): DataFrame containing the observations.
        date_col (str): The column name with the observation dates in ISO
            format.
        gap_days (int): Number of days to define a gap. Observations with gaps
            greater than this will start a new group.

    Returns:
        pd.DataFrame: DataFrame with an additional 'group' column indicating
        sequential groups.
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


def parse_output_filename(lon, lat, day_of_year, observtn, obs_group):
    """
    Generate a filename for sky background data based on coordinates, day of
    the year, and observation details.

    Args:
        lon (int or float): Longitude (0 to 360).
        lat (int or float): Latitude (-90 to 90).
        day_of_year (int): Day of the year (1 to 365).
        observtn (int): Observation number.
        obs_group (int): Observation group number.

    Returns:
        str: Formatted filename string for the sky background FITS file.
    """
    try:
        # Ensure inputs are integers where appropriate
        lon = int(lon)
        lat = int(lat)
        day_of_year = int(day_of_year)
        observtn = int(observtn)
        obs_group = int(obs_group)
    except ValueError:
        raise ValueError("All inputs must be convertible to integers.")

    # Construct the sky position encoding
    sky_encoder = (
        f"LON{lon:03d}_LAT{lat:+03d}_DOY{day_of_year:03d}_"
        f"{observtn:03}_{obs_group:03d}"
    )

    return f"nis_clear_gr700xd_{sky_encoder}_skyflat.fits"


def write_fits_file(keys, bkg, data, err, masks, metadata_table, output_dir):

    bkg_denoise = apply_denoise_nl_means(bkg)

    # pad images to re-stitch reference pixels
    padwidth = (4, 4)
    bkg = np.pad(bkg, (padwidth, padwidth))
    bkg_denoise = np.pad(bkg_denoise, (padwidth, padwidth))
    data = np.pad(data, ((0, 0), padwidth, padwidth))
    masks = np.pad(
        masks,
        ((0, 0), padwidth, padwidth),
        mode="constant",
        constant_values=True,
    )

    program, targprop, observtn, obs_group = keys

    n_samps = len(data)

    day_of_year = metadata_table["day_of_year"].mean()
    ecliptic_lon = metadata_table["ecliptic_lon"].mean()
    ecliptic_lat = metadata_table["ecliptic_lat"].mean()

    outfile = parse_output_filename(
        ecliptic_lon, ecliptic_lat, day_of_year, observtn, obs_group
    )

    primary = fits.PrimaryHDU()
    primary.header["DATE"] = get_current_utc(), "Date file created"
    primary.header["FILENAME"] = outfile
    primary.header["AUTHOR"] = "Tyler Baines"
    primary.header["DESCRIP"] = "Empirical SOSS Background"
    primary.header["FLATCORR"] = True, "Flatfield corrected images"
    primary.header["NSAMPS"] = (n_samps, "# of images for Emp. background")
    primary.header["PROGRAM"] = program, "Program ID"
    primary.header["TARGPROP"] = targprop, "Target Property Id."
    primary.header["LON"] = ecliptic_lon, "Ecliptic longitude"
    primary.header["LAT"] = ecliptic_lat, "Ecliptic latitude"
    primary.header["DAY"] = np.median(day_of_year), "median day of year"
    primary.header["MINDAY"] = np.min(day_of_year), "min day of year"
    primary.header["MAXDAY"] = np.max(day_of_year), "max day of year"

    # Background ImageHDU
    bkg_hdu = fits.ImageHDU(data=bkg, name="BKG")
    bkg_hdu.header["DESCRIP"] = "Flatfielded Empirical SOSS Background"

    # Smooth image using Non-local means denoising algorithm.

    bkg_dn_hdu = fits.ImageHDU(data=bkg_denoise, name="BKGDN")
    bkg_dn_hdu.header["DESCRIP"] = (
        "Flatfielded Empirical SOSS Background NLM denoised"
    )

    # ImageHDU for masked array of the input images w/ source masks
    image_hdu = fits.ImageHDU(
        data=np.ma.array(data=data_list, mask=source_masks),
        name="IMAGES",
    )

    outfits = fits.HDUList([primary, bkg_hdu, bkg_dn_hdu, image_hdu])

    outfits.writeto(output_dir / outfile, overwrite=True)


if __name__ == "__main__":
    # environment variables
    zodical_csv = os.getenv("ZODICAL_CSV")
    flatfield_dir = Path(os.getenv("FLATFIELD_DIR"))
    sossbkg_dir = Path(os.getenv("SOSSBKG_DIR"))

    parallel = True
    overwrite = False

    print(f"Flatfield Dir: {flatfield_dir}")
    print(f"SOSSBkg Dir: {sossbkg_dir}")

    # load data
    zodical_stage1_df = pd.read_csv(zodical_csv)
    zodical_stage1_df = zodical_stage1_df.sort_values(["PROGRAM", "DATE-OBS"])

    # Ensure the output directory exists
    flatfield_dir.mkdir(parents=True, exist_ok=True)

    # Update the DataFrame with the processed file paths
    fits_files = zodical_stage1_df["FILEPATH"]

    # apply flatfield step to all files
    apply_flatfield(
        files=fits_files,
        output_dir=flatfield_dir,
        parallel=parallel,
        overwrite=overwrite,
    )

    zodical_stage1_df["FLATFIELD_FILEPATH"] = zodical_stage1_df[
        "FILEPATH"
    ].apply(
        lambda x: flatfield_dir
        / Path(x).name.replace(".fits", "_flatfield.fits")
    )
    # Step for running through files to group them to make bkg images

    # keys which to group data by
    group_keys = ["PROGRAM", "TARGPROP", "OBSERVTN"]

    include_pwcpos = False

    if include_pwcpos:
        group_keys += ["PWCPOS"]

    sossbkg_dir.mkdir(exist_ok=True, parents=True)

    # group sequential observations that are under a day
    zodical_stage1_df["obs_groups"] = zodical_stage1_df.groupby(group_keys)[
        "DATE-BEG"
    ].transform(group_sequential_observations)

    # grouped dataframes
    grouped_dfs = zodical_stage1_df.groupby(group_keys + ["obs_groups"])

    # Main loop to iterate over the grouped data frames by the Program ID,
    # target, and observation number. Datasets are then given a group ID based
    # on observations that were taking sequentially. There are a total of 10
    # exposuresper observation for the Calibration acitivty programs that
    # should have been taken one after another. However, im guessing due to
    # scheduling windows some data was taken a few or many days apart. We know
    # that the background intensity changes over the course of the year. In
    # addition, it's hypothesized that the background may vary spatailly given
    # the PWCPOS keyword.

    for keys, df in grouped_dfs:

        # iterate over the observation groups to process and combine the
        # file datasets
        file_list = df["FLATFIELD_FILEPATH"]  # .to_list()

        # Load the data into JWST DataModel
        input_models = load_fits_data(file_list)

        # metadata extraction
        metadata_table = parse_metadata_table(input_models)
        metadata_table["FILEPATH"] = fits_files

        # Extract data, error, and dq arrays from input models
        data_list = np.array([dm.data[0][nonref] for dm in input_models])
        err_list = np.array([dm.err[0][nonref] for dm in input_models])
        dq_list = np.array([dm.dq[0][nonref] for dm in input_models])

        # intial estimate of the background median combined sigma-clipped
        _, bkg, _ = sigma_clipped_stats(data_list, sigma=1, axis=0)

        #  Source Masking step
        source_masks = create_source_masks(
            data_list, err_list, dq_list, bkg, nsigma=2
        )

        # Estimate the background again using apply the source masks
        _, bkg, _ = sigma_clipped_stats(
            data_list,
            mask=source_masks,
            sigma=1,
            axis=0,
        )

        # Fill NaNs
        bkg = apply_maskfill(bkg, mask=np.isnan(bkg))

        # Lastly, write results to a file
        write_fits_file(
            bkg, data_list, err_list, source_masks, metadata_table, sossbkg_dir
        )
