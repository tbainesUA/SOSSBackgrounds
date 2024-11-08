import multiprocessing

from pathlib import Path

from jwst.flatfield import FlatFieldStep


def process_flatfield(file, output_dir):
    """
    Applies the flat field correction to a given file and saves the result.

    Parameters:
    - file (str or Path): Path to the input FITS file.
    - output_dir (str or Path): Directory where the corrected file will be
    saved.

    Returns:
    - Path: The output file path where the corrected FITS file is saved.
    """
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)

    try:
        flat_field_step = FlatFieldStep()

        # Process the file using the FlatFieldStep
        result = flat_field_step.call(file)

        # Construct the new filename and output path
        new_filename = result.meta.filename.replace(".fits", "_flatfield.fits")
        outfile = output_dir / new_filename

        # Save the result to the output directory
        result.save(outfile)

        return outfile

    except Exception as e:
        print(f"Error processing file {file}: {e}")
        return None


# Apply the flat field correction to each file in the DataFrame
# Adding error handling and progress tracking
def apply_flat_field(files, output_dir, parallel=False):
    """
    Applies the flat field correction to a DataFrame containing file paths.

    Parameters:
    - df (pd.DataFrame): DataFrame containing file paths.
    - filepath_col (str): Column name in the DataFrame that contains the file 
        paths.
    - output_dir (Path): Directory where processed files will be saved.

    Returns:
    - pd.Series: Series containing the paths to the processed files.
    """

    output_dirs = [output_dir for _ in range(len(files))]

    if parallel:
        ncores = multiprocessing.cpu_count() // 2
        with multiprocessing.Pool(ncores) as pool:
            processed_filepaths = pool.starmap(
                process_flatfield, zip(files, output_dirs)
            )
    else:
        processed_filepaths = []
        for i, file in enumerate(files):
            print(f"Processing file {i+1}/{len(files)}: {file}")
            processed_file = process_flatfield(file, output_dir)
            processed_filepaths.append(processed_file)

    return processed_filepaths
