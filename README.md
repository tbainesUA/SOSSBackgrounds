

# SOSSBackgrounds

## Characterization of the GR700XD Sky Background for NIRISS/SOSS Observations

This repository provides tools for creating empirical sky background models for NIRISS/SOSS observations. The code is set up to utilize internal STScI network storage. For external users, datasets can be queried and acquired using astroquery.

## Overview
This repository contains scripts to generate empirical SOSS backgrounds using data from Programs 4479 and 6658. These are dedicated calibration programs that collect SOSS background observations in FULL mode for a variety of sky positions and multiple dither pointings. In addition to these datasets, other FULL mode observations are included.

The generated backgrounds will be used to improve the background subtraction step in the SOSS data reduction process. Current background subtraction techniques, which rely on a static background template, often result in residuals of a few percent when scaling the template across different observations. This method does not perform optimally, and our work aims to refine this process.

### Key Points:
* Program 4479 & 6658: These programs provide the calibration data for SOSS background observations.
* FULL Mode: Observations are conducted in FULL mode across various sky positions and dither pointings.
* Goal: Improve background subtraction for SOSS data reduction by creating more accurate, empirical background models.