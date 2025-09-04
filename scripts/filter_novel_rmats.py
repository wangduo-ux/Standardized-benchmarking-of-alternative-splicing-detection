import os
import pandas as pd

def filter_rmats_novel_events(df, software, event, input_dir, sample_name):
    """
    Filter out novel or de novo alternative splicing events based on the software type.

    Parameters:
    - df: pd.DataFrame, input AS event data.
    - software: str, software name ('rMATS', 'MAJIQ', 'PSI-Sigma', 'SUPPA2').
    - event: str, event type ('SE', 'A3SS', 'A5SS', 'RI', 'MX').
    - input_dir: str, base input directory containing software-specific folders.
    - sample_name: str, sample identifier.

    Returns:
    - pd.DataFrame: filtered AS events with novel/de novo events removed.
    """
    if "rmats" in software.lower():
        # Define key columns to identify unique events for each event type
        event_key_columns = {
            "SE": ['GeneID','chr', 'strand', 'exonStart_0base', 'exonEnd', 
                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'],
            "A3SS": ['GeneID','chr', 'strand', 'longExonStart_0base', 'longExonEnd', 
                     'flankingES', 'flankingEE', 'shortES', 'shortEE'],
            "A5SS": ['GeneID','chr', 'strand', 'longExonStart_0base', 'longExonEnd', 
                     'flankingES', 'flankingEE', 'shortES', 'shortEE'],
            "RI": ['GeneID','chr', 'strand', 'riExonStart_0base', 'riExonEnd', 
                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE'],
            "MX": ['GeneID','chr', 'strand', '1stExonStart_0base', '1stExonEnd', 
                   '2ndExonStart_0base', '2ndExonEnd', 'upstreamES', 'upstreamEE', 
                   'downstreamES', 'downstreamEE']
        }

        key_columns = event_key_columns.get(event)
        if key_columns is None:
            print(f"Unsupported event type for rMATS: {event}")
            return df

        # Construct file paths for novel splice site and junction events
        folder = os.path.join(input_dir, software, sample_name)
        splice_file = os.path.join(folder, f"fromGTF.novelSpliceSite.{event}.txt")
        junction_file = os.path.join(folder, f"fromGTF.novelJunction.{event}.txt")

        # Load files if they exist
        dfs_to_remove = []
        for f in [splice_file, junction_file]:
            if os.path.exists(f):
                try:
                    df_tmp = pd.read_csv(f, sep='\t', dtype=str)
                    dfs_to_remove.append(df_tmp)
                except Exception as e:
                    print(f"Failed to load {f}: {e}")

        if not dfs_to_remove:
            return df

        # Combine and merge with input df to filter out novel events
        df_remove = pd.concat(dfs_to_remove, ignore_index=True)
        common_keys = [col for col in key_columns if col in df.columns and col in df_remove.columns]
        if not common_keys:
            print(f"[Warning] {event} has no matching columns with novel events")
            return df

        # Ensure matching columns are string type for merging
        for col in common_keys:
            df[col] = df[col].astype(str)
            df_remove[col] = df_remove[col].astype(str)

        df_merged = df.merge(df_remove[common_keys].drop_duplicates(), 
                             on=common_keys, how='left', indicator=True)
        df_filtered = df_merged[df_merged['_merge'] == 'left_only'].drop(columns=['_merge'])

    elif "majiq" in software.lower():
        # Remove events labeled as denovo
        event_ids_to_remove = df.loc[df['denovo'] == True, 'event_id'].unique()
        df_filtered = df[~df['event_id'].isin(event_ids_to_remove)]

    elif "psi-sigma" in software.lower():
        # Keep only known events
        df_filtered = df[df['type'] != "novel"]

    elif "suppa2" in software.lower():
        # SUPPA2 does not annotate novel events separately; keep all
        df_filtered = df

    return df_filtered
