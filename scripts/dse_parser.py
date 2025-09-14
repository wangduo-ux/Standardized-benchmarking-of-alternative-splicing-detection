import pandas as pd
import numpy as np
import os
import glob
from typing import Tuple, Set

def annotate_dse_class(df: pd.DataFrame, event_type: str, software: str) -> pd.DataFrame:
    """
    Annotate DSE class (up-regulate, down-regulate, non-DSE) based on software-specific criteria.
    """
    df = df.copy()
    
    if 'rmats' in software.lower():
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
        df['IncLevelDifference'] = pd.to_numeric(df['IncLevelDifference'], errors='coerce')
        conditions = [
            (df['FDR'] <= 0.01) & (df['IncLevelDifference'] >= 0.05),
            (df['FDR'] <= 0.01) & (df['IncLevelDifference'] <= -0.05)
        ]
        choices = ['up-regulate', 'down-regulate']
        df['class'] = np.select(conditions, choices, default='non-DSE')
    
    elif 'suppa2' in software.lower():
        dpsi_col = [col for col in df.columns if 'dPSI' in col][0]
        pval_col = [col for col in df.columns if 'p-val' in col][0]
        conditions = [
            (df[dpsi_col] >= 0.05) & (df[pval_col] < 0.05),
            (df[dpsi_col] <= -0.05) & (df[pval_col] < 0.05)
        ]
        choices = ['up-regulate', 'down-regulate']
        df['class'] = np.select(conditions, choices, default='non-DSE')
        
    elif 'majiq' in software.lower():
        prob_cols = [col for col in df.columns if 'probability_changing' in col]
        dpsi_cols = [col for col in df.columns if 'median_dpsi' in col]
        df = df.dropna(subset=prob_cols + dpsi_cols, how='any')
        for uid, group in df.groupby('uniform_ID'):
            if event_type == "SE":
                group_inclusion = group[group['spliced_with'] == 'A']
            elif event_type in ["A3SS", "A5SS"]:
                group_inclusion = group[group['junction_name'] == 'Proximal']
            elif event_type in ["AF", "AL"]:
                group_inclusion = group[group['junction_name'] == 'Distal']
            elif event_type == "MX":
                strand = group['strand'].iloc[0]
                group_inclusion = group[group['junction_name'].isin(
                    ["C1_A1", "C2_A2"] if strand == '+' else ["C1_A2", "C2_A1"]
                )]
            elif event_type == "RI":
                group_inclusion = group[group['junction_name'].str.contains('intron')]
            else:
                continue

            if group_inclusion.empty:
                df.loc[df['uniform_ID'] == uid, 'class'] = 'non-DSE'
                continue

            prob_high = (group_inclusion[prob_cols] > 0.95).any(axis=None)
            max_prob_row = group_inclusion.loc[group_inclusion[prob_cols].max(axis=1).idxmax()]
            dpsi_val = float(max_prob_row[dpsi_cols[0]])

            if prob_high:
                df.loc[df['uniform_ID'] == uid, 'class'] = 'up-regulate' if dpsi_val > 0 else 'down-regulate'
            else:
                df.loc[df['uniform_ID'] == uid, 'class'] = 'non-DSE'

    elif 'psi-sigma' in software.lower():
        if event_type == "MX":
            df = df.loc[df.groupby('uniform_ID')['Target Exon'].idxmin()]
        df['ΔPSI (%)'] = pd.to_numeric(df['ΔPSI (%)'], errors='coerce')
        df['FDR (BH)'] = pd.to_numeric(df['FDR (BH)'], errors='coerce')
        conditions = [
            (df['ΔPSI (%)'] >= 5) & (df['FDR (BH)'] < 0.05),
            (df['ΔPSI (%)'] <= -5) & (df['FDR (BH)'] < 0.05)
        ]
        choices = ['up-regulate', 'down-regulate']
        df['class'] = np.select(conditions, choices, default='non-DSE')

    else:
        raise ValueError(f"Unsupported software: {software}")

    return df


def ase_extract(df: pd.DataFrame, software: str, event_type: str, sample_name: str, input_dir: str) -> Tuple[Set[str], Set[str]]:
    """
    Extract ASE (allele-specific expression) IDs for test and control groups.

    Returns
    -------
    Tuple of sets: (test_ase_ids, control_ase_ids)
    """
    psivec = df.copy()

    if 'rmats' in software.lower():
        psivec['mean_psi_test'] = psivec['IncLevel1'].apply(lambda x: np.mean([float(i) for i in str(x).split(',') if i != 'NA']) if pd.notnull(x) else np.nan)
        psivec['mean_psi_control'] = psivec['IncLevel2'].apply(lambda x: np.mean([float(i) for i in str(x).split(',') if i != 'NA']) if pd.notnull(x) else np.nan)
        test_ase_ids = psivec.loc[(psivec['mean_psi_test'] >= 0.05) & (psivec['mean_psi_test'] <= 0.95), 'uniform_ID'].tolist()
        control_ase_ids = psivec.loc[(psivec['mean_psi_control'] >= 0.05) & (psivec['mean_psi_control'] <= 0.95), 'uniform_ID'].tolist()

    elif 'suppa2' in software.lower():
        event_map = {"SE": "SE", "A3SS": "A3", "A5SS": "A5", "AF": "AF", "AL": "AL", "RI": "RI", "MX": "MX"}
        event = event_map.get(event_type)
        search_path = os.path.join(input_dir, 'SUPPA2', sample_name, '*psivec')
        matched_files = glob.glob(search_path)
        psivec = pd.read_csv(matched_files[0], sep="\t")
        psivec = psivec[psivec.index.str.contains(event)]
        columns = psivec.columns.tolist()
        sample_groups = list(set(['_'.join(col.split('_')[:2]) for col in columns]))
        psivec['mean_psi_test'] = psivec.loc[:, psivec.columns.str.contains(sample_groups[0])].mean(axis=1)
        psivec['mean_psi_control'] = psivec.loc[:, psivec.columns.str.contains(sample_groups[1])].mean(axis=1)
        test_ase_ids = psivec.loc[(psivec['mean_psi_test'] >= 0.05) & (psivec['mean_psi_test'] <= 0.95)].index.tolist()
        control_ase_ids = psivec.loc[(psivec['mean_psi_control'] >= 0.05) & (psivec['mean_psi_control'] <= 0.95)].index.tolist()

    elif 'majiq' in software.lower():
        # similar filtering as in annotate_dse_class
        if event_type == "SE":
            psivec = psivec[psivec['spliced_with'] == 'A']
        elif event_type in ["A3SS", "A5SS"]:
            psivec = psivec[psivec['junction_name'] == 'Proximal']
        elif event_type in ["AF", "AL"]:
            psivec = psivec[psivec['junction_name'] == 'Distal']
        elif event_type == "MX":
            psivec = psivec[((psivec['strand'] == '+') & (psivec['junction_name'].isin(["C1_A1", "C2_A2"]))) |
                            ((psivec['strand'] == '-') & (psivec['junction_name'].isin(['C1_A2', "C2_A1"])))]
        elif event_type == "RI":
            psivec = psivec[psivec['junction_name'].str.contains('intron')]
        psi_cols = [col for col in psivec.columns if 'median_psi' in col]
        test_ase_ids = (psivec.groupby('uniform_ID')[psi_cols[1]]
                        .apply(lambda group: group.between(0.05, 0.95).any().any())
                        .loc[lambda x: x].index.tolist())
        control_ase_ids = (psivec.groupby('uniform_ID')[psi_cols[0]]
                           .apply(lambda group: group.between(0.05, 0.95).any().any())
                           .loc[lambda x: x].index.tolist())

    elif 'psi-sigma' in software.lower():
        psivec['mean_psi_test'] = psivec['T Values'].apply(lambda x: np.mean([float(i) for i in str(x).split('|') if i.lower() != 'na']) if pd.notnull(x) else np.nan)
        psivec['mean_psi_control'] = psivec['N Values'].apply(lambda x: np.mean([float(i) for i in str(x).split('|') if i.lower() != 'na']) if pd.notnull(x) else np.nan)
        test_ase_ids = (psivec.groupby('uniform_ID')['mean_psi_test']
                        .apply(lambda group: group.between(5, 95).any().any())
                        .loc[lambda x: x].index.tolist())
        control_ase_ids = (psivec.groupby('uniform_ID')['mean_psi_control']
                           .apply(lambda group: group.between(5, 95).any().any())
                           .loc[lambda x: x].index.tolist())

    else:
        raise ValueError(f"Unsupported software: {software}")

    return set(test_ase_ids), set(control_ase_ids)
