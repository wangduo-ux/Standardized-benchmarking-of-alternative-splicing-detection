import pandas as pd
import numpy as np


def annotate_dse_class(df: pd.DataFrame, event_type, software: str) -> pd.DataFrame:
    """
    Annotate differential splicing events (DSE) based on
    software-specific statistical thresholds.

    Returns a DataFrame with columns:
        - uniform_ID
        - class (DSE / non-DSE)
    """

    df = df.copy()
    software_lower = software.lower()

    # ===== rMATS =====
    if 'rmats' in software_lower:
        df['FDR'] = pd.to_numeric(df['FDR'], errors='coerce')
        df['IncLevelDifference'] = pd.to_numeric(df['IncLevelDifference'], errors='coerce')

        df['class'] = 'non-DSE'
        df.loc[
            (df['FDR'] <= 0.01) &
            (df['IncLevelDifference'].abs() >= 0.05),
            'class'
        ] = 'DSE'

    # ===== SUPPA2 =====
    elif 'suppa2' in software_lower:
        dpsi_col = [col for col in df.columns if 'dPSI' in col][0]
        pval_col = [col for col in df.columns if 'p-val' in col][0]

        df['class'] = 'non-DSE'
        df.loc[
            (df[pval_col] <= 0.05) &
            (df[dpsi_col].abs() >= 0.05),
            'class'
        ] = 'DSE'

    # ===== MAJIQ =====
    elif 'majiq' in software_lower:
        df['class'] = 'non-DSE'
        df.loc[df["event_changing"], 'class'] = 'DSE'

    # ===== PSI-Sigma =====
    elif 'psi-sigma' in software_lower:
        df['ΔPSI (%)'] = pd.to_numeric(df['ΔPSI (%)'], errors='coerce')
        df['FDR (BH)'] = pd.to_numeric(df['FDR (BH)'], errors='coerce')

        cond = (
            (df['ΔPSI (%)'].abs() >= 10) &
            (df['FDR (BH)'] < 0.01)
        )

        if event_type == "MX":
            df['class'] = np.where(
                cond.groupby(df['uniform_ID']).transform('any'),
                'DSE',
                'non-DSE'
            )
        else:
            df['class'] = np.where(cond, 'DSE', 'non-DSE')

    # ===== Spladder =====
    elif 'spladder' in software_lower:
        df['p_val_adj'] = pd.to_numeric(df['p_val_adj'], errors='coerce')
        df['dPSI'] = pd.to_numeric(df['dPSI'], errors='coerce')

        df['class'] = 'non-DSE'
        df.loc[
            (df['p_val_adj'] <= 0.05) &
            (df['dPSI'].abs() >= 0.05),
            'class'
        ] = 'DSE'

    # ===== Whippet =====
    elif 'whippet' in software_lower:
        df['Probability'] = pd.to_numeric(df['Probability'], errors='coerce')
        df['DeltaPsi'] = pd.to_numeric(df['DeltaPsi'], errors='coerce')

        df['class'] = 'non-DSE'
        df.loc[
            (df['Probability'] >= 0.9) &
            (df['DeltaPsi'].abs() >= 0.01),
            'class'
        ] = 'DSE'

    else:
        raise ValueError(f"Unsupported software: {software}")

    # ===== final output =====
    df = df[["uniform_ID", "class"]].drop_duplicates()

    return df
