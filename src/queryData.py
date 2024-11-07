import os
from pathlib import Path
import pandas as pd
from dotenv import load_dotenv


def load_csv(file_path: str) -> pd.DataFrame:
    """Loads a CSV file into a DataFrame."""
    return pd.read_csv(file_path)


def filter_dataframe(df: pd.DataFrame, query_criteria: str) -> pd.DataFrame:
    """Filters a DataFrame based on query criteria."""
    return df.query(query_criteria)


def save_dataframes(dataframes: dict):
    """Saves DataFrames to CSV files in the specified directory."""

    for filename, df in dataframes.items():
        df.to_csv(filename, index=False)


if __name__ == "__main__":

    # Load environment variables
    load_dotenv()

    CSV_DIR = Path(os.getenv("CSV_DIR"))

    CSV_DIR.mkdir(parents=True, exist_ok=True)

    # internal soss library database
    SOSSLIB_STAGE1_CSV = os.getenv("SOSSLIB_STAGE1_CSV")
    SOSSLIB_UNCAL_CSV = os.getenv("SOSSLIB_UNCAL_CSV")

    # output file
    ZODICAL_CSV = os.getenv("ZODICAL_CSV")
    ZODICAL_VALIDATION_CSV = os.getenv("ZODICAL_VALIDATION_CSV")
    ZODICAL_SCI_CSV = os.getenv("ZODICAL_SCI_CSV")

    # Load SOSS library data files
    uncal_df = load_csv(SOSSLIB_UNCAL_CSV)
    stage1_df = load_csv(SOSSLIB_STAGE1_CSV)

    # Define common filtering criteria
    subarrays = ["FULL"]
    filters = ["CLEAR"]
    zodiCal_programs = [1541, 4479, 6658]

    # query criteria to filter soss observations
    query = (
        "SUBARRAY == @subarrays and "
        "FILTER == @filters and "
        "PROGRAM >= 1091 and "
        "(FILENAME.str.contains('uncal') or FILENAME.str.contains('rateints'))"
    )

    # Filter data
    uncal_df_filtered = filter_dataframe(uncal_df, query)
    stage1_df_filtered = filter_dataframe(stage1_df, query)

    # Create DataFrames for different purposes
    zodical_stage1_df = stage1_df_filtered.query("PROGRAM in @zodiCal_programs")
    zodical_stage1_validation_df = stage1_df_filtered.query(
        "PROGRAM not in @zodiCal_programs"
    )
    zodical_uncal_sci_df = uncal_df_filtered.query("CATEGORY in ['GTO', 'GO']")

    # Define output data and files
    output_files = {
        ZODICAL_CSV: zodical_stage1_df,
        ZODICAL_VALIDATION_CSV: zodical_stage1_validation_df,
        ZODICAL_SCI_CSV: zodical_uncal_sci_df,
    }

    # Save filtered data to CSV files
    save_dataframes(output_files)
