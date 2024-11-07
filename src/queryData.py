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


def save_dataframes(output_dir: Path, dataframes: dict):
    """Saves DataFrames to CSV files in the specified directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, df in dataframes.items():
        df.to_csv(output_dir / filename, index=False)


if __name__ == "__main__":

    # Load environment variables
    load_dotenv()
    
    CSV_DIR = Path(os.getenv("CSV_DIR"))

    SOSSLIB_STAGE1_CSV = os.getenv("SOSSLIB_STAGE1_CSV")
    SOSSLIB_UNCAL_CSV = os.getenv("SOSSLIB_UNCAL_CSV")

    # Load SOSS library data files
    uncal_df = load_csv(SOSSLIB_UNCAL_CSV)
    stage1_df = load_csv(SOSSLIB_STAGE1_CSV)

    # Define common filtering criteria
    subarrays = ["FULL"]
    filters = ["CLEAR"]
    zodiCal_programs = [1541, 4479, 6658]

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
        "zodical_stage1_rateints.csv": zodical_stage1_df,
        "zodical_stage1_validation_rateints.csv": zodical_stage1_validation_df,
        "zodical_uncal_sci.csv": zodical_uncal_sci_df,
    }

    # Save filtered data to CSV files
    save_dataframes(CSV_DIR, output_files)