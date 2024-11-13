import multiprocessing
from pathlib import Path
from jwst.flatfield import FlatFieldStep


def process_flatfield(file, output_dir, overwrite=False):
    """
    Applies the flat field correction to a given file and saves the result.

    Parameters:
    - file (str or Path): Path to the input FITS file.
    - output_dir (str or Path): Directory where the corrected file will be saved.
    - overwrite (bool): If True, reprocess the file even if it already exists.

    Returns:
    - Path or None: The output file path where the corrected FITS file is saved,
                    or None if the file was skipped.
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)

    # Construct the output filename and check if it exists
    new_filename = Path(file).stem + "_flatfield.fits"
    outfile = output_dir / new_filename

    # Check if the file already exists before any processing is done
    if outfile.exists() and not overwrite:
        print(f"File {outfile} already exists. Skipping (overwrite=False).")
        return outfile  # Skip processing if file exists and overwrite is False

    try:
        # Process the file using the FlatFieldStep only if needed
        flat_field_step = FlatFieldStep()
        result = flat_field_step.process(file)

        # Save the result to the output directory
        result.save(outfile)
        print(f"Processed and saved: {outfile}")
        return outfile

    except Exception as e:
        print(f"Error processing file {file}: {e}")
        return None


def apply_flatfield(files, output_dir, parallel=False, overwrite=False):
    """
    Applies the flat field correction to a list of files.

    Parameters:
    - files (list): List containing paths to FITS files.
    - output_dir (Path): Directory where processed files will be saved.
    - parallel (bool): If True, use multiprocessing to process files in parallel.
    - overwrite (bool): If True, reprocess files even if they already exist.

    Returns:
    - list: List containing the paths to the processed files.
    """
    output_dirs = [output_dir for _ in range(len(files))]
    overwrite_flags = [overwrite for _ in range(len(files))]

    if parallel:
        ncores = multiprocessing.cpu_count() // 2
        with multiprocessing.Pool(ncores) as pool:
            processed_filepaths = pool.starmap(
                process_flatfield, zip(files, output_dirs, overwrite_flags)
            )
    else:
        processed_filepaths = []
        for i, file in enumerate(files):
            print(f"Processing file {i+1}/{len(files)}: {file}")
            processed_file = process_flatfield(file, output_dir, overwrite)
            processed_filepaths.append(processed_file)

    return processed_filepaths
