import pandas as pd
import numpy as np

def select_rows(group):
    mask = group['class'].isin(['up-regulate','down-regulate'])
    if mask.any():
        return group[mask] 
    else:
        return group.iloc[[0]]  
      
def annotate_dse(df: pd.DataFrame,event_type, software: str) -> pd.DataFrame:
    """
    Annotate DSE class (up-regulate, down-regulate, non-DSE) based on software-specific criteria.
    """
    df = df.copy()
    
    if 'rmats'in software.lower():
        # Ensure numeric conversion
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
        pval_col = [col for col in df.columns if 'p-val' in col in col][0]
        conditions = [
        (df[dpsi_col] >= 0.05) & (df[pval_col] < 0.05),
        (df[dpsi_col] <= -0.05) & (df[pval_col] < 0.05)
        ]
        choices = ['up-regulate', 'down-regulate']
        df['class'] = np.select(conditions, choices, default='non-DSE')
        
    elif  'majiq' in software.lower():
        prob_cols = [col for col in df.columns if 'probability_changing' in col]
        dpsi_cols = [col for col in df.columns if 'median_dpsi' in col]
        dpsi_col = dpsi_cols[0]   # 'lab4_8-lab4_6_median_dpsi'
        prob_col = prob_cols[0]
        df = df.dropna(subset=[prob_col, dpsi_col], how="any")
        df = df.copy()
        if event_type == "SE":
            df = df[df['spliced_with'] == 'A']
            conditions = [(df[dpsi_col] > 0) & (df[prob_col] > 0.95),(df[dpsi_col] < 0) & (df[prob_col] > 0.95)]
            choices = ["up-regulate", "down-regulate"]
            df["class"] = np.select(conditions, choices, default="non-DSE")
            df = df.groupby('uniform_ID', group_keys=False).apply(select_rows)
              
        elif event_type in ["A3SS", "A5SS"]:
            df = df[df['junction_name'] == 'Proximal']
            conditions = [(df[dpsi_col] > 0) & (df[prob_col] > 0.95),(df[dpsi_col] < 0) & (df[prob_col] > 0.95)]
            choices = ["up-regulate", "down-regulate"]
            df["class"] = np.select(conditions, choices, default="non-DSE")

        elif event_type in ["AF", "AL"]:
            df = df[df['junction_name'] == 'Distal']
            conditions = [(df[dpsi_col] > 0) & (df[prob_col] > 0.95),(df[dpsi_col] < 0) & (df[prob_col] > 0.95)]
            choices = ["up-regulate", "down-regulate"]
            df["class"] = np.select(conditions, choices, default="non-DSE")
            
        elif event_type == "MX":
            df_positive = df[df["strand"]=="+"]  
            df_positive = df_positive[df_positive['junction_name'].isin(["C1_A1", "C2_A2"])]
            df_negative = df[df["strand"]=="-"]
            df_negative = df_negative[df_negative['junction_name'].isin(["C1_A2", "C2_A1"])]
            df = pd.concat([df_positive, df_negative], ignore_index=True)
            conditions = [(df[dpsi_col] > 0) & (df[prob_col] > 0.95),(df[dpsi_col] < 0) & (df[prob_col] > 0.95)]
            choices = ["up-regulate", "down-regulate"]
            df["class"] = np.select(conditions, choices, default="non-DSE")
            df = df.groupby('uniform_ID', group_keys=False).apply(select_rows)

        elif event_type == "RI":
            df = df[df['junction_name'].str.contains('intron')]
            conditions = [(df[dpsi_col] > 0) & (df[prob_col] > 0.95),(df[dpsi_col] < 0) & (df[prob_col] > 0.95)]
            choices = ["up-regulate", "down-regulate"]
            df["class"] = np.select(conditions, choices, default="non-DSE")

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
    df = df[["uniform_ID","class"]]
    return df  
