import os
import glob
import pandas as pd

from .gtf_parser import extract_gene_transcript_strand

def load_file(input_dir, sample_name, event_type, software, groupA, groupB, gtf_file):
    """
    Load splicing results from different software tools and return a standardized DataFrame.
    """

    df = None  # ensure defined

    # ================= SUPPA2 =================
    if "suppa2" in software.lower():
        event_type_map = {
            "SE": "SE", "A3SS": "A3", "A5SS": "A5",
            "AF": "AF", "AL": "AL", "RI": "RI", "MX": "MX"
        }

        type_in_suppa2 = event_type_map.get(event_type)

        search_path = os.path.join(input_dir, sample_name,software, "*dpsi")
        matched_files = glob.glob(search_path)

        if len(matched_files) == 1:
            df = pd.read_csv(matched_files[0], sep="\t", comment="#")
            df["event_id"] = df.index

            df = df[
                df["event_id"].astype(str)
                .str.split(";", n=1).str[1]
                .str.contains(type_in_suppa2, na=False)
            ].copy()

        elif len(matched_files) > 1:
            raise FileNotFoundError(f"Multiple dpsi files found: {search_path}")
        else:
            raise FileNotFoundError(f"No dpsi file found: {search_path}")

    # ================= rMATS =================
    elif "rmats" in software.lower():
        event = "MXE" if event_type == "MX" else event_type

        file_path = os.path.join(input_dir, sample_name,software, f"{event}.MATS.JCEC.txt")

        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        df = pd.read_csv(file_path, sep="\t", comment="#")

    # ================= psi-sigma =================
    elif "psi-sigma" in software.lower():
        search_path = os.path.join(input_dir, sample_name,software, "*sorted.txt")
        matched_files = glob.glob(search_path)
        # Load GTF annotation
        gtf_df = extract_gene_transcript_strand(gtf_file)
        #search_path_2 = os.path.join(input_dir, software, sample_name, "*mapping.txt")
        #matched_files_2 = glob.glob(search_path_2)

        if len(matched_files) != 1:
            raise FileNotFoundError(f"Invalid sorted.txt file: {search_path}")

        #if len(matched_files_2) != 1:
            #raise FileNotFoundError(f"Invalid mapping.txt file: {search_path_2}")

        df = pd.read_csv(matched_files[0], sep="\t", comment="#")

        df = df.rename(columns={"Reference Transcript": "transcript_id"})

        df["type"] = df["transcript_id"].astype(str).str.split(".").str[0].apply(
            lambda x: "novel" if x == "Ex" else "-"
        )

        df["transcript_id"] = (
            df["transcript_id"]
            .astype(str)
            .str.replace(r"^((Ex|TSS)\.)+", "", regex=True)
        )

        #anno_df = pd.read_csv(matched_files_2[0], sep="\t", comment="#")
        #anno_df.columns = ["transcript_id", "gene_id", "strand"]

        if event_type == "SE":
            df = df[df["Event Type"] == "SES"].copy()
        elif event_type == "A3SS":
            df = df[df["Event Type"].isin(["A3SS", "TSS|A3SS"])].copy()
        elif event_type == "A5SS":
            df = df[df["Event Type"].isin(["A5SS", "TSS|A5SS"])].copy()
        elif event_type == "MX":
            df = df[df["Event Type"] == "MXS"].copy()
        elif event_type == "RI":
            df = df[df["Event Type"].isin(["IR", "IR (overlapping region)"])].copy()
        
        df["key"] = df["transcript_id"].str.replace(r"[^A-Za-z0-9]", "", regex=True)
        gtf_df["key"] = gtf_df["transcript_id"].str.replace(r"[^A-Za-z0-9]", "", regex=True)
        df = pd.merge(df, gtf_df, on="key", how="left")

    # ================= MAJIQ =================
    elif "majiq" in software.lower():
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
            raise ValueError(f"Unsupported MAJIQ event: {event_type}")

        if isinstance(files, str):
            df = pd.read_csv(os.path.join(input_dir, sample_name,software, files),
                             sep="\t", comment="#", low_memory=False)
        else:
            dfs = [
                pd.read_csv(os.path.join(input_dir, sample_name,software, f),
                            sep="\t", comment="#", low_memory=False)
                for f in files
            ]
            df = pd.concat(dfs, axis=0)

    # ================= SplAdder =================
    elif "spladder" in software.lower():
        event_map = {
            "SE": "exon_skip",
            "A3SS": "alt_3prime",
            "A5SS": "alt_5prime",
            "RI": "intron_retention",
            "MX": "mutex_exons",
        }

        if event_type not in event_map:
            raise ValueError(f"Unsupported event_type for SplAdder: {event_type}")

        prefix = event_map[event_type]
        base_dir = os.path.join(input_dir, sample_name)

        tsv_pattern = os.path.join(base_dir,software,"*testing*",f"*extended*C[0-9]_{prefix}*.tsv")
        tsv_files = glob.glob(tsv_pattern)

        if not tsv_files:
            raise FileNotFoundError(f"SplAdder result not found: {tsv_pattern}")

        df = pd.read_csv(tsv_files[0], sep="\t", comment="#").iloc[:, :15].copy()

        merge_pattern = os.path.join(base_dir,software, f"*graphs_{prefix}_C[0-9].confirmed.txt.gz")
        merge_files = glob.glob(merge_pattern)

        if not merge_files:
            raise FileNotFoundError(f"Merge graph file not found: {prefix}")

        merge_graphs = pd.read_csv(merge_files[0], sep="\t", comment="#", low_memory=False)

        psi_cols_groupA = [c for c in merge_graphs.columns if "psi" in c.lower() and any(g in c for g in groupA)]
        psi_cols_groupB = [c for c in merge_graphs.columns if "psi" in c.lower() and any(g in c for g in groupB)]

        renameA = {c: f"test_psi_{i}" for i, c in enumerate(psi_cols_groupA, 1)}
        renameB = {c: f"control_psi_{i}" for i, c in enumerate(psi_cols_groupB, 1)}

        merge_graphs.rename(columns={**renameA, **renameB}, inplace=True)

        df_subset = merge_graphs[
            ["strand", "event_id", "is_annotated", *renameA.values(), *renameB.values()]
        ]

        df = pd.merge(df, df_subset, on="event_id", how="left")

    # ================= Whippet =================
    elif "whippet" in software.lower():
        search_path = os.path.join(input_dir, sample_name,software, '*diff.gz')
        matched_files = glob.glob(search_path)
        if len(matched_files) != 1:
            raise FileNotFoundError(f"no or more than one diff.gz file in the {search_path}")
        file_path = matched_files[0]
        df = pd.read_csv(file_path, sep="\t", comment="#", index_col=False)
        if event_type == "SE":
            df = df[df['Type']=='CE'].copy()
        elif event_type == "A3SS":    
            df = df[df['Type']=='AA'].copy()
        elif event_type == "A5SS":  
            df = df[df['Type']=='AD'].copy()
        elif event_type == "AF":
            df = df[df['Type']=='AF'].copy()
        elif event_type == "AL":
            df = df[df['Type']=='AL'].copy()
        elif event_type == "RI":
            df = df[df['Type']=='RI'].copy()

    # ================= fallback =================
    else:
        raise ValueError(f"Unsupported software: {software}")

    return df
  
  
def save_dataframe(df, output_dir, sample_name, software, event_type):
    """
    Save a DataFrame to a standardized output file.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to be saved.
    output_dir : str
        Root output directory.
    sample_name : str
        Sample name.
    software : str
        Software/tool name.
    event_type : str
        Splicing event type.
    """

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Create subdirectory: output_dir/software/sample_name/
    sub_dir = os.path.join(output_dir, software, sample_name)
    os.makedirs(sub_dir, exist_ok=True)

    # Output file name
    output_file = os.path.join(
        sub_dir,
        f"{software}.{event_type}.uid.txt"
    )

    # Save file
    df.to_csv(output_file, sep="\t", index=False)
