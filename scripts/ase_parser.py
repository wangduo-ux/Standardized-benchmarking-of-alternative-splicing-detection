import os
import glob
import numpy as np
import pandas as pd

from .parser import PARSER_FUNCTIONS


def ase_extract(df, software, event_type, sample_name, input_dir):
    """
    Extract ASE (Alternative Splicing Events) for test and control groups
    from different software outputs.
    """

    psivec = df.copy()

    if 'rmats' in software.lower():
        psivec['mean_psi_test'] = psivec['IncLevel1'].apply(
            lambda x: np.mean([float(i) for i in str(x).split(',') if i != 'NA'])
            if pd.notnull(x) else np.nan
        )
        psivec['mean_psi_control'] = psivec['IncLevel2'].apply(
            lambda x: np.mean([float(i) for i in str(x).split(',') if i != 'NA'])
            if pd.notnull(x) else np.nan
        )

        test_ase_df = psivec.loc[
            (psivec['mean_psi_test'] >= 0.05) & (psivec['mean_psi_test'] <= 0.95),
            ["uniform_ID", "mean_psi_test"]
        ]

        control_ase_df = psivec.loc[
            (psivec['mean_psi_control'] >= 0.05) & (psivec['mean_psi_control'] <= 0.95),
            ["uniform_ID", "mean_psi_control"]
        ]

    elif 'suppa2' in software.lower():

        event_map = {
            "SE": "SE", "A3SS": "A3", "A5SS": "A5",
            "AF": "AF", "AL": "AL", "RI": "RI", "MX": "MX"
        }

        type_in_suppa2 = event_map.get(event_type)

        search_path = os.path.join(input_dir,  sample_name,software, '*psivec')
        matched_files = glob.glob(search_path)

        psivec = pd.read_csv(matched_files[0], sep="\t")
        columns = psivec.columns.tolist()

        psivec["event_id"] = psivec.index
        psivec = psivec[
            psivec["event_id"].str.split(";", n=1).str[1].str.contains(type_in_suppa2, na=False)
        ].copy()

        parser = PARSER_FUNCTIONS["SUPPA2"][event_type]
        psivec = parser(psivec, event_type)

        half = len(columns) // 2
        test_cols = columns[:half]
        ctrl_cols = columns[half:]

        psivec["mean_psi_test"] = psivec[test_cols].mean(axis=1)
        psivec["mean_psi_control"] = psivec[ctrl_cols].mean(axis=1)

        test_ase_df = psivec.loc[
            (psivec['mean_psi_test'] >= 0.05) & (psivec['mean_psi_test'] <= 0.95),
            ["uniform_ID", "mean_psi_test"]
        ]

        control_ase_df = psivec.loc[
            (psivec['mean_psi_control'] >= 0.05) & (psivec['mean_psi_control'] <= 0.95),
            ["uniform_ID", "mean_psi_control"]
        ]

    elif 'majiq' in software.lower():

        psi_cols = [col for col in psivec.columns if 'median_psi' in col]
        test_psi_col = psi_cols[1]
        control_psi_col = psi_cols[0]

        test_ase_df = (
            psivec
            .groupby("uniform_ID")[test_psi_col]
            .apply(lambda x: ",".join(map(str, x)) if ((x >= 0.05) & (x <= 0.95)).any() else None)
            .dropna()
            .reset_index()
            .rename(columns={test_psi_col: "mean_psi_test"})
        )

        control_ase_df = (
            psivec
            .groupby("uniform_ID")[control_psi_col]
            .apply(lambda x: ",".join(map(str, x)) if ((x >= 0.05) & (x <= 0.95)).any() else None)
            .dropna()
            .reset_index()
            .rename(columns={control_psi_col: "mean_psi_control"})
        )

    elif 'psi-sigma' in software.lower():

        psivec['mean_psi_test'] = psivec['T Values'].apply(
            lambda x: np.mean([float(i) for i in str(x).split('|') if i.lower() != 'na'])
            if pd.notnull(x) else np.nan
        )

        psivec['mean_psi_control'] = psivec['N Values'].apply(
            lambda x: np.mean([float(i) for i in str(x).split('|') if i.lower() != 'na'])
            if pd.notnull(x) else np.nan
        )

        control_ase_df = (
            psivec
            .groupby("uniform_ID")["mean_psi_control"]
            .apply(lambda x: ",".join(map(str, x)) if x.between(5, 95).any() else None)
            .dropna()
            .reset_index(name="mean_psi_control")
        )

        test_ase_df = (
            psivec
            .groupby("uniform_ID")["mean_psi_test"]
            .apply(lambda x: ",".join(map(str, x)) if x.between(5, 95).any() else None)
            .dropna()
            .reset_index(name="mean_psi_test")
        )

    elif 'spladder' in software.lower():

        psivec["mean_psi_test"] = psivec.filter(like="test_psi").mean(axis=1)
        psivec["mean_psi_control"] = psivec.filter(like="control_psi").mean(axis=1)

        test_ase_df = psivec.loc[
            (psivec['mean_psi_test'] >= 0.05) & (psivec['mean_psi_test'] <= 0.95),
            ["uniform_ID", "mean_psi_test"]
        ]

        control_ase_df = psivec.loc[
            (psivec['mean_psi_control'] >= 0.05) & (psivec['mean_psi_control'] <= 0.95),
            ["uniform_ID", "mean_psi_control"]
        ]

    elif 'whippet' in software.lower():

        test_ase_df = (
            psivec.loc[
                (psivec['Psi_A'] >= 0.05) & (psivec['Psi_A'] <= 0.95),
                ["uniform_ID", "Psi_A"]
            ]
            .rename(columns={"Psi_A": "mean_psi_test"})
        )

        control_ase_df = (
            psivec.loc[
                (psivec['Psi_B'] >= 0.05) & (psivec['Psi_B'] <= 0.95),
                ["uniform_ID", "Psi_B"]
            ]
            .rename(columns={"Psi_B": "mean_psi_control"})
        )

    else:
        raise ValueError(f"Unsupported software: {software}")

    return test_ase_df, control_ase_df
