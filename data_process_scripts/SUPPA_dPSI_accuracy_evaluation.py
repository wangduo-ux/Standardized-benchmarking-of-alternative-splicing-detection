#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Alternative splicing (AS) evaluation script based on SUPPA results

This script evaluates:
1. Reproducibility across technical replicates (CV-like metric)
2. Accuracy of dPSI estimation against truth
3. Optional inter-lab concordance (disabled by default)

All results are printed to stdout.
"""

import argparse
import os
import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression


# =========================
# Argument parsing
# =========================
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lab", required=True, help="Lab prefix, e.g. lab17_")
    parser.add_argument("--suppa_dir", required=True, help="SUPPA output base directory")
    parser.add_argument("--truth", required=True, help="AS truth file")
    parser.add_argument("--truth_pcr", required=True, help="PCR truth file")
    return parser.parse_args()


# =========================
# Load truth
# =========================
def load_truth(truth_file, truth_pcr_file):
    truth = pd.read_csv(truth_file)
    truth["ASE_compare"] = truth["compare"] + "_" + truth["ASE"]

    truth_pcr = pd.read_table(truth_pcr_file)
    truth_pcr["ASE_compare"] = truth_pcr["sample_pair"] + "_" + truth_pcr["ASE"]
    truth_pcr = truth_pcr.rename(columns={"delta_psi": "mean_delta_psi_mean"})

    return truth, truth_pcr


# =========================
# Load PSI replicates
# =========================
def load_psi(lab,suppa_dir):
    dfs = []
    for sample in ["M8", "F7", "D5", "D6"]:
        psi_file = f"{suppa_dir}/{lab}{sample}.psi"
        if not os.path.exists(psi_file):
            continue

        df = pd.read_table(psi_file)
        df = df.reset_index()
        df.columns = ["ASE", "psi_1", "psi_2", "psi_3"]
        df["ASE_compare"] = sample + "_" + df["ASE"]
        dfs.append(df)

    df = pd.concat(dfs, ignore_index=True).dropna()
    df["mean"] = df[["psi_1", "psi_2", "psi_3"]].mean(axis=1)

    # remove near-constitutive events
    df = df[abs(df["mean"] - 0.5) <= 0.45]
    return df


# =========================
# Reproducibility metric
# =========================
def cv_metric(df):
    u = df[["psi_1", "psi_2", "psi_3"]].mean(axis=1)
    s = df[["psi_1", "psi_2", "psi_3"]].std(axis=1, ddof=0)
    G = len(df)
    return np.sqrt(s.sum() / G)


# =========================
# Load dPSI
# =========================
def load_dpsi(lab, suppa_dir):
    dfs = []
    for sample in ["M8", "F7", "D5"]:
        f = f"{suppa_dir}/{lab}{sample}.dpsi"
        if not os.path.exists(f):
            continue

        df = pd.read_table(f).reset_index()
        df.columns = ["ASE", "dpsi", "p"]
        df["ASE_compare"] = sample + "/D6_" + df["ASE"]
        dfs.append(df[["ASE_compare", "dpsi"]])

    return pd.concat(dfs, ignore_index=True).dropna()


# =========================
# Accuracy metrics
# =========================
def accuracy(df):
    pearson, _ = scipy.stats.pearsonr(df["dpsi"], df["mean_delta_psi_mean"])
    spearman, _ = spearmanr(df["dpsi"], df["mean_delta_psi_mean"])
    rmse = np.sqrt(((df["dpsi"] - df["mean_delta_psi_mean"]) ** 2).mean())
    #nrmse = rmse / (df["dpsi"].max() - df["dpsi"].min())

    X = df[["mean_delta_psi_mean"]]
    y = df["dpsi"]

    return pearson, spearman, rmse


# =========================
# Main
# =========================
def main():
    args = parse_args()

    truth, truth_pcr = load_truth(args.truth, args.truth_pcr)

    psi = load_psi(args.lab, args.suppa_dir)
    im = cv_metric(psi)

    dpsi = load_dpsi(args.lab, args.suppa_dir)
    merged = pd.merge(dpsi, truth, on="ASE_compare", how="inner")
    merged_pcr = pd.merge(dpsi, truth_pcr, on="ASE_compare", how="inner")

    metrics = accuracy(merged)
    metrics_pcr = accuracy(merged_pcr)
    print(
        "\t".join(
            map(
                str,
                [
                    args.lab,
                    round(im, 4),
                    len(merged),
                    *[round(x, 4) for x in metrics],
                    len(merged),
                    *[round(x, 4) for x in metrics_pcr],
                ],
            )
        )
    )


if __name__ == "__main__":
    main()
