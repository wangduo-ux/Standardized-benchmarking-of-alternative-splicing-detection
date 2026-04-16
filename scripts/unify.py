import os
import re
import glob
import itertools
from itertools import combinations
import subprocess 
import numpy as np
import pandas as pd

# ===== import your internal functions =====
# =========================
from .io import (load_file,save_dataframe)


from .ase_parser import ase_extract

from .dse_annotation import annotate_dse_class

from .dpsi_parser import extract_dpsi

from .parser import (
    parse_software
)

# pairwise + event integration
from .pairwise_analysis import (
    compute_pairwise_pearson,
    summarize_whippet_overlap,
    build_pairwise_df,
    build_upset_df
)
from .event_merger import merge_splicing_events

def get_support_software_list(event_type):
    if event_type in ('SE', 'A3SS', 'A5SS', 'RI'):
        return ['SUPPA2', 'rMATS', 'PSI-Sigma', 'MAJIQ','Spladder','Whippet']
    elif event_type in ('AF', 'AL'):
        return ['SUPPA2', 'MAJIQ','Whippet']
    elif event_type in ('MX'):
        return ['SUPPA2', 'rMATS', 'PSI-Sigma', 'MAJIQ', 'Spladder'] 
    else:
        return []
      
def unify_results(input_dir, software_list, event_list, sample_name,
                  output_dir, gtf_file, groupA, groupB):
    """
    Iterate over each software and event type, parse matching files,
    and concatenate all results.
    """

    # ===== Create output directories =====
    os.makedirs(output_dir, exist_ok=True)
    tmp_folder = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_folder, exist_ok=True)

    # ===== Locate R scripts (relative path) =====
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    R_DIR = os.path.join(BASE_DIR, "scripts", "r_scripts")

    dse_event_dict = {}
    test_ase_dict = {}
    control_ase_dict = {}
    dse_gene_dict = {}
    summary_all = []
    dpsi_list = []

    for ev in event_list:
        if ev not in ["SE", "A3SS", "A5SS", "MX", "RI", "AF", "AL"]:
            print(f"Warning: unsupported event type '{ev}', skipping.")
            continue

        support_software = get_support_software_list(ev)

        query_software = [
            sw2 for sw2 in software_list
            if any(sw1.lower() in sw2.lower() for sw1 in support_software)
        ]

        test_ase_dfs = {}
        control_ase_dfs = {}
        dse_dfs = {}
        dpsi_dfs = {}

        for sw in query_software:
            print(f"Running event type: {ev}, tool: {sw}")

            df = load_file(input_dir, sample_name, ev, sw, groupA, groupB, gtf_file)

            df = parse_software(df, sw, ev)
            
            if df is None or df.empty:
                continue

            # ===== Save parsed result =====
            save_dataframe(df, output_dir, sample_name, sw, ev)

            # ===== ASE extraction =====
  
            test_ase_dfs[sw], control_ase_dfs[sw] = ase_extract(
                df, sw, ev, sample_name, input_dir
            )

            # ===== DSE annotation =====
            dse_df = annotate_dse_class(df, ev, sw)
            dse_dfs[sw] = dse_df[dse_df["class"] == "DSE"]

            num_DSE = dse_df[dse_df["class"] == "DSE"]["uniform_ID"].nunique()
            num_nonDSE = dse_df[dse_df["class"] == "non-DSE"]["uniform_ID"].nunique()

            summary_all.append({
                "Software": sw,
                "Event_type": ev,
                "DSE": num_DSE,
                "nonDSE": num_nonDSE
            })

            # ===== dPSI extraction =====
            dpsi_dfs[sw] = extract_dpsi(df, sw, ev)

        # ---------- ASE merge ----------
        test_ase = merge_splicing_events(test_ase_dfs, ev)
        upset_df = build_upset_df(test_ase)
        upset_df["event_type"] = ev
        test_ase_dict[ev] = upset_df

        control_ase = merge_splicing_events(control_ase_dfs, ev)
        upset_df = build_upset_df(control_ase)
        upset_df["event_type"] = ev
        control_ase_dict[ev] = upset_df

        # ---------- DSE merge ----------
        dse = merge_splicing_events(dse_dfs, ev)
        upset_df = build_upset_df(dse)
        upset_df["event_type"] = ev
        dse_event_dict[ev] = upset_df

        # ---------- dPSI merge ----------
        dpsi = merge_splicing_events(dpsi_dfs, ev)
        dpsi = dpsi[["event_id", "tool", "uniform_ID", "value", "type"]]
        dpsi_list.append(dpsi)

        # ---------- DSG ----------
        dsg = dse[["tool", "type", "gene"]]
        dse_gene_dict[ev] = dsg

    # ===== Final outputs =====

    # ASE (test)
    test_ase_df = pd.concat(test_ase_dict.values(), ignore_index=True)
    test_ase_file = os.path.join(tmp_folder, f"{sample_name}_test_ase_upset_plot.csv")
    test_ase_df.to_csv(test_ase_file, sep="\t", index=False)

    whippet_test = summarize_whippet_overlap(test_ase_df)
    if isinstance(whippet_test, pd.DataFrame):
        whippet_test_file = os.path.join(tmp_folder, f"{sample_name}_whippet_test_ase_overlap_plot.csv")
        whippet_test.to_csv(whippet_test_file, sep="\t", index=False)

    # ASE (control)
    control_ase_df = pd.concat(control_ase_dict.values(), ignore_index=True)
    control_file = os.path.join(tmp_folder, f"{sample_name}_control_ase_upset_plot.csv")
    control_ase_df.to_csv(control_file, sep="\t", index=False)

    whippet_control = summarize_whippet_overlap(control_ase_df)
    if isinstance(whippet_control, pd.DataFrame):
        whippet_control_file = os.path.join(tmp_folder, f"{sample_name}_whippet_control_ase_overlap_plot.csv")
        whippet_control.to_csv(whippet_control_file, sep="\t", index=False)

        subprocess.run([
            "Rscript",
            os.path.join(R_DIR, "ASB_plot4.R"),
            whippet_control_file,
            os.path.join(output_dir, f"{sample_name}_whippet_overlap.pdf")
        ], check=True)

    # DSE
    dse_df = pd.concat(dse_event_dict.values(), ignore_index=True)
    dse_file = os.path.join(tmp_folder, f"{sample_name}_dse_upset_plot.csv")
    dse_df.to_csv(dse_file, sep="\t", index=False)

    # Summary
    summary_df = pd.DataFrame(summary_all)
    summary_file = os.path.join(tmp_folder, f"{sample_name}_dse_number_summary.csv")
    summary_df.to_csv(summary_file, sep="\t", index=False)

    # dPSI
    dpsi_df = pd.concat(dpsi_list, ignore_index=True)
    final = build_pairwise_df(dpsi_df)
    dpsi_file = os.path.join(tmp_folder, f"{sample_name}_dpsi_plot.csv")
    final.to_csv(dpsi_file, sep="\t", index=False)

    # ===== R plotting =====
    subprocess.run(["Rscript", os.path.join(R_DIR, "ASB_plot1.R"), test_ase_file,
                    os.path.join(output_dir, f"{sample_name}_test_ase_upsetplot.pdf")], check=True)

    subprocess.run(["Rscript", os.path.join(R_DIR, "ASB_plot1.R"), control_file,
                    os.path.join(output_dir, f"{sample_name}_control_ase_upsetplot.pdf")], check=True)

    subprocess.run(["Rscript", os.path.join(R_DIR, "ASB_plot1.R"), dse_file,
                    os.path.join(output_dir, f"{sample_name}_dse_upsetplot.pdf")], check=True)

    subprocess.run(["Rscript", os.path.join(R_DIR, "ASB_plot2.R"), summary_file,
                    os.path.join(output_dir, f"{sample_name}_dse_num_plot.pdf")], check=True)

    subprocess.run(["Rscript", os.path.join(R_DIR, "ASB_plot3.R"), dpsi_file,
                    os.path.join(output_dir, f"{sample_name}_dpsi.pdf")], check=True)

    compute_pairwise_pearson(final, output_dir, sample_name)
