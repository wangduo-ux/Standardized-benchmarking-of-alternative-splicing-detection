# utils.py (or file_io.py)
import os
import pandas as pd

def save_dataframe(df: pd.DataFrame, output_dir: str, sample_name: str, software: str, event_type: str) -> None:
    """
    Save a DataFrame to a tab-delimited file in a structured directory.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to save.
    output_dir : str
        Root directory to save output files.
    sample_name : str
        Name of the sample.
    software : str
        Name of the software (e.g., SUPPA2, rMATS, MAJIQ, PSI-Sigma).
    event_type : str
        Type of splicing event (e.g., SE, A3SS, A5SS, MX, RI, AF, AL).

    Returns
    -------
    None
    """
    if not os.path.exists(output_dir):
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")

    sub_dir = os.path.join(output_dir, software, sample_name)
    os.makedirs(sub_dir, exist_ok=True)

    output_file = os.path.join(sub_dir, f"{software}.{event_type}.uid.txt")
    df.to_csv(output_file, sep='\t', index=False)
