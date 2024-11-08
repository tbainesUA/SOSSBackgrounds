import numpy as np
from scipy import ndimage

from photutils.segmentation import detect_sources, detect_threshold


def detect_source_mask(data, err=None, mask=None, bkg=None, nsigma=2):
    """Source Detection step"""
    im = data.copy()

    # Mask data
    if mask is None:
        mask = np.isnan(data)

    # Detection step
    thresh = detect_threshold(
        data=im, nsigma=nsigma, background=bkg, error=err, mask=mask
    )
    seg = detect_sources(im, thresh, npixels=10, connectivity=8)

    # apply dialation of any found pixels as neighbors could be affected too.
    source_mask = ndimage.binary_dilation(seg.data > 0)

    # combine with the input mask
    source_mask = np.logical_or(mask, source_mask)

    return source_mask
