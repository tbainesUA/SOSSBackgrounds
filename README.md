

# SOSSBackgrounds

## Characterization of the GR700XD Sky Background for NIRISS/SOSS Observations

This repository provides tools for creating empirical sky background models for NIRISS/SOSS observations. The code is set up to utilize internal STScI network storage. For external users, datasets can be queried and acquired using astroquery.

## Overview
This repository contains scripts to generate empirical SOSS backgrounds using data from Programs 4479 and 6658. These are dedicated calibration programs that collect SOSS background observations in FULL mode for a variety of sky positions and multiple dither pointings. In addition to these datasets, other FULL mode observations are included.

The generated backgrounds will be used to improve the background subtraction step in the SOSS data reduction process. Current background subtraction techniques, which rely on a static background template, often result in residuals of a few percent when scaling the template across different observations. This method does not perform optimally, and our work aims to refine this process.

Output filename: nis_clear_gr700xd_LON009_LAT-02_DOY320_001_000_skyflat.fits

### Key Points:
* Program 4479 & 6658: These programs provide the calibration data for SOSS background observations.
* FULL Mode: Observations are conducted in FULL mode across various sky positions and dither pointings.
* Goal: Improve background subtraction for SOSS data reduction by creating more accurate, empirical background models.


---

### Output Files and File Naming Scheme

The scripts in this repository generate background model files in FITS format. These files represent the sky background observations and are used in the SOSS data reduction process. Each output file corresponds to a specific sky position and observation, which are reflected in the file name.

All output files include the following extensions:

  1.  BKG    : Empirical Sci. Background Data
  2.  BKGDN  : Empirical Sci. Background Data (Smoothed/Denoised)
  3.  IMAGES(rename to SCI) : Images used to Produce Empirical Background   
  4.  MASKS  : Source Masks
  
to include in outputs:

  5.  ERR : Error Array 
  6.  DQ  : Data Quality Array 

### File Naming Scheme

The output files follow a structured naming convention to encode essential metadata, including sky position, day of the year, observation number, and group number. The naming scheme is designed to be concise while providing relevant information about each file.

#### Format:

```
nis_clear_gr700xd_LON{lon:03d}_LAT{lat:+03d}_DOY{day_of_year:03d}_{observtn:03}_{obs_group:03d}_skyflat.fits
```
#### Where:

* nis: Denotes the NIRISS detector.
* clear: Refers to the CLEAR filter used in the observations.
* gr700xd: Specifies the grism used for the observation.
* LON{lon:03d}: Encodes the ecliptic longitude (0–360) of the sky position, formatted as a 3-digit integer.
* LAT{lat:+03d}: Encodes the ecliptic latitude (-90 to 90), formatted as a signed 3-digit integer (e.g., LAT+45, LAT-10).
* DOY{day_of_year:03d}: The day of the year (1–365), formatted as a 3-digit integer.
* {observtn:03}: The observation number, formatted as a 3-digit integer.
* {obs_group:03d}: The observation group number, formatted as a 3-digit integer.
* skyflat.fits: A suffix indicating that the file contains a sky background model (also referred to as a "skyflat").

#### Example:

An observation at ecliptic longitude 90° and latitude +45°. Taken on the 123rd day of the year (DOY123). The 45th observation in the 1st group. 
```
nis_clear_gr700xd_LON090_LAT+45_DOY123_045_001_skyflat.fits
```


#### Purpose of Output Files

The generated `.fits` files represent sky background models derived from the specified observations. These models will be used for background subtraction during the SOSS data reduction process, ensuring more accurate results across different observations and reducing residual errors.