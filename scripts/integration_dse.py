import os
import pandas as pd
from functools import reduce
from scripts.dse_annotation import annotate_dse

def load_files(input_dir, event_list, software, sample_name):
    """
    Load AS event files for the specified software and samples, returning merged DSE dataframe.
    """
    software_lower = software.lower()
    if 'suppa2' in software_lower or 'majiq' in software_lower:
        support_event = list(set(event_list) & set(["SE", "A3SS", "A5SS", "AF", "AL", "RI", "MX"]))
    elif 'rmats' in software_lower or 'psi-sigma' in software_lower:
        support_event = list(set(event_list) & set(["SE", "A3SS", "A5SS", "RI", "MX"]))
    else:
        raise ValueError(f"Unsupported software: {software}")
        
    dse_list = []
    
    for ev in support_event:
        file_path = os.path.join(input_dir, software, sample_name, f"{software}.{ev}.uid.txt")
        if not os.path.exists(file_path):
            print(f"Warning: file not found: {file_path}")
            continue
            
        temp_df = pd.read_csv(file_path, sep="\t", low_memory=False)
        dse = annotate_dse(temp_df, ev, software)
        dse_list.append(dse)
            
    dse_df = pd.concat(dse_list, axis=0, ignore_index=True)
    return dse_df

def integration_analysis(input_dir, software_list, event_list, sample_name, output_dir):
    """
    Integrate multiple software results and output summary tables.
    """
    from functools import reduce
    import os
    
    dfs = []
    for sw in software_list:
        dse_df = load_files(input_dir, event_list, sw, sample_name)
        dse_df['class'] = dse_df['class'].replace(['up-regulate', 'down-regulate'], 'DSE')
        dse_df.rename(columns={'class': f'{sw}'}, inplace=True)
        dfs.append(dse_df)
    
    df_combined = reduce(lambda left, right: pd.merge(left, right, on='uniform_ID', how='outer'), dfs)
    dse_cols = df_combined.columns[1:1+len(software_list)]
    df_combined["sum"] = df_combined[dse_cols].isin(["DSE"]).sum(axis=1)
    
    num_software = len(software_list)
    for n in range(num_software):
        df_combined['class'] = ['DSE' if x > n else 'non-DSE' for x in df_combined['sum']]
        sub_dir = os.path.join(output_dir, sample_name) 
        os.makedirs(sub_dir, exist_ok=True)
        support_software = n + 1
        output_file = os.path.join(sub_dir, f"integration_{support_software}of{n}.txt")
        df_combined.to_csv(output_file, sep='\t', index=False)
