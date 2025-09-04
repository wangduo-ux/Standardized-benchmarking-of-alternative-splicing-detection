# file_loader.py
import os
import glob
import pandas as pd
from typing import Union

def load_file(input_dir: str, sample_name: str, event_type: str, software: str) -> pd.DataFrame:
    """
    Load splicing event file for a given software, sample, and event type.

    Parameters:
    -----------
    input_dir : str
        Root directory containing software output folders.
    sample_name : str
        Name of the sample to load.
    event_type : str
        Splicing event type, e.g., 'SE', 'A3SS', 'A5SS', 'MX', 'RI', 'AF', 'AL'.
    software : str
        Software name: 'SUPPA2', 'rMATS', 'MAJIQ', or 'PSI-Sigma'.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the splicing event data.
    """
    software_lower = software.lower()

    if "suppa2" in software_lower:
        search_path = os.path.join(input_dir, 'SUPPA2', sample_name, '*dpsi*')
        matched_files = glob.glob(search_path)
        if not matched_files:
            raise FileNotFoundError(f"No dpsi file found in {search_path}")
        file_path = matched_files[0]
        df = pd.read_csv(file_path, sep='\t', comment='#')

    elif "rmats" in software_lower:
        event = "MXE" if event_type == "MX" else event_type
        file_path = os.path.join(input_dir, 'rMATS', sample_name, f'{event}.MATS.JCEC.txt')
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        df = pd.read_csv(file_path, sep='\t', comment='#')

    elif "psi-sigma" in software_lower:
        file_path = os.path.join(input_dir, 'PSI-Sigma', sample_name, 'PSIsigma_r10_ir3.sorted.txt')
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        df = pd.read_csv(file_path, sep="\t", comment='#')

    elif "majiq" in software_lower:
        majiq_file_map = {
            "SE": "cassette.tsv",
            "A3SS": ("alt3prime.tsv", "p_alt3prime.tsv"),
            "A5SS": ("alt5prime.tsv", "p_alt5prime.tsv"),
            "AF": ("alternate_first_exon.tsv", "p_alternate_first_exon.tsv"),
            "AL": ("alternate_last_exon.tsv", "p_alternate_last_exon.tsv"),
            "RI": "alternative_intron.tsv",
            "MX": "mutually_exclusive.tsv"
        }
        files = majiq_file_map.get(event_type)
        if files is None:
            raise ValueError(f"No files mapped for event type: {event_type}")

        if isinstance(files, str):
            df = pd.read_csv(os.path.join(input_dir, "MAJIQ", sample_name, files),
                             sep="\t", comment="#", low_memory=False)
        else:
            dfs = [pd.read_csv(os.path.join(input_dir, "MAJIQ", sample_name, f),
                               sep="\t", comment="#", low_memory=False) for f in files]
            df = pd.concat(dfs, axis=0)

    else:
        raise ValueError(f"Unknown software: {software}")

    return df
