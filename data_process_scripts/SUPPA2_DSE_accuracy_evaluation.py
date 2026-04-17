#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SUPPA DAS accuracy evaluation (single lab)

This script evaluates:
1. DAS classification accuracy (TP / FP / FN / TN)
2. Precision / Recall / F1 / MCC
3. DAS recovery & inconsistency between technical replicates

All results are printed to stdout.
"""

import argparse
import math
import pandas as pd
import numpy as np


# ======================
# Argument parsing
# ======================
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--suppa_dir", required=True,
                        help="SUPPA result directory")
    parser.add_argument("--lab", required=True,
                        help="Lab ID prefix, e.g. lab1_")
    parser.add_argument("--truth", required=True,
                        help="Ground truth DAS annotation file")
    parser.add_argument("--truth_pcr", required=True,
                        help="PCR validated DAS truth file")
    return parser.parse_args()


# ======================
# DAS annotation rule
# ======================
def annotate_das(row):
    if row["p"] < 0.05 and row["dpsi"] >= 0.05:
        return "up-regulate"
    elif row["p"] < 0.05 and row["dpsi"] <= -0.05:
        return "down-regulate"
    else:
        return "non-DAS"


# ======================
# Load dPSI results
# ======================
def load_dpsi(lab, suppa_dir, sample):
    file_path = f"{suppa_dir}/{lab}{sample}.dpsi"
    df = pd.read_table(file_path).reset_index()

    df.columns = ["ASE", "dpsi", "p"]
    df["ASE"] = f"{sample}/D6_" + df["ASE"]
    df["anno"] = df.apply(annotate_das, axis=1)

    return df


# ======================
# Classification
# ======================
def classify(df):
    conditions = [
        (df["final"] == "up-regulate") & (df["anno"] == "up-regulate"),
        (df["final"] == "up-regulate") & (df["anno"] != "up-regulate"),

        (df["final"] == "down-regulate") & (df["anno"] == "down-regulate"),
        (df["final"] == "down-regulate") & (df["anno"] != "down-regulate"),

        (df["final"] == "non-DAS") & (df["anno"] == "non-DAS"),
        (df["final"] == "non-DAS") & (df["anno"] != "non-DAS"),
    ]

    values = ["TP", "FN", "TP", "FN", "TN", "FP"]
    df["class"] = np.select(conditions, values, default="FP")

    return df


# ======================
# Metrics calculation
# ======================
def calculate_metrics(df):
    TP = (df["class"] == "TP").sum()
    FP = (df["class"] == "FP").sum()
    FN = (df["class"] == "FN").sum()
    TN = (df["class"] == "TN").sum()

    recall = TP / (TP + FN) if TP + FN else 0
    precision = TP / (TP + FP) if TP + FP else 0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall else 0
    )

    denom = math.sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
    mcc = (TP * TN - FP * FN) / denom if denom else 0

    return TP, FP, TN, FN, recall, precision, f1, mcc


# ======================
# Main
# ======================
def main():
    args = parse_args()

    # Load truth annotations
    truth_pcr = pd.read_table(args.truth_pcr)
    truth_pcr["ASE"] = truth_pcr["sample_pair"] + "_" + truth_pcr["ASE"]
    truth_pcr = truth_pcr.rename(columns={"DAS": "final"})
    
    truth = pd.read_csv(args.truth)
    truth["ASE"] = truth["compare"] + "_" + truth["ASE"]
    truth = truth.rename(columns={"DAS": "final"})

    # Load SUPPA results
    dpsi_all = pd.concat([
        load_dpsi(args.lab,  args.suppa_dir, s)
        for s in ["M8", "F7", "D5"]
    ])

    # Merge with truth
    merged = pd.merge(dpsi_all, truth, on="ASE", how="inner")
    merged = classify(merged)
    
    merged_pcr = pd.merge(dpsi_all, truth_pcr, on="ASE", how="inner")
    merged_pcr = classify(merged_pcr)
    # Metrics
    TP, FP, TN, FN, recall, precision, f1, mcc = calculate_metrics(merged)
    TP_2, FP_2, TN_2, FN_2, recall_2, precision_2, f1_2, mcc_2 = calculate_metrics(merged_pcr)
    # Output
    print(
    "lab\t"
    "Quartet_truth\tTP\tFP\tTN\tFN\trecall\tprecision\tF1\tMCC\t"
    "RTqPCR_truth\tTP\tFP\tTN\tFN\trecall\tprecision\tF1\tMCC"
    )
    print(
        "\t".join(
            map(
                str,
                [
                  args.lab,"Quartet reference",TP, FP, TN, FN, recall, precision, f1, mcc,"RT-qPCR reference",TP_2, FP_2, TN_2, FN_2, recall_2, precision_2, f1_2, mcc_2],
            )
        )
    )


if __name__ == "__main__":
    main()
