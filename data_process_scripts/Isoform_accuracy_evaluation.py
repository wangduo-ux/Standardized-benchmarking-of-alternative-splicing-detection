#!/usr/bin/env python

import numpy as np
import pandas as pd
import scipy.stats
import argparse


# --------------------------------------------------
# Core processing functions
# --------------------------------------------------
def process(df):
    """
    Calculate log2FC from 6 columns:
    first 3 vs last 3
    """
    df = df.reset_index(drop=True)

    # log2 transform
    df.iloc[:, 1:7] = np.log2(df.iloc[:, 1:7] + 0.01)

    df["mean1"] = df.iloc[:, 1:4].mean(axis=1)
    df["mean2"] = df.iloc[:, 4:7].mean(axis=1)

    df["log2FC"] = df["mean1"] - df["mean2"]

    return df.loc[:, ["transcript_id", "log2FC"]]


def processdata(df, lab):
    """
    Generate three comparisons and concatenate
    """

    data1 = df[
        ["transcript_id",
         f"{lab}202207", f"{lab}202208", f"{lab}202209",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    data2 = df[
        ["transcript_id",
         f"{lab}202210", f"{lab}202211", f"{lab}202212",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    data3 = df[
        ["transcript_id",
         f"{lab}202213", f"{lab}202214", f"{lab}202215",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    data1 = process(data1)
    data2 = process(data2)
    data3 = process(data3)

    data1["transcript_id"] = "M8/D6_" + data1["transcript_id"].astype(str)
    data2["transcript_id"] = "F7/D6_" + data2["transcript_id"].astype(str)
    data3["transcript_id"] = "D5/D6_" + data3["transcript_id"].astype(str)

    return pd.concat([data1, data2, data3], axis=0)


# --------------------------------------------------
# Evaluation
# --------------------------------------------------
def evaluate_fc_accuracy(df_subset, lab, truth, truth2):
    """
    Evaluate FC accuracy using RMSE and Pearson correlation
    """

    processed_df = processdata(df_subset, lab)

    combine1 = pd.merge(truth, processed_df, on="transcript_id", how="inner")
    combine2 = pd.merge(truth2, processed_df, on="transcript_id", how="inner")

    combine2 = combine2[
        ~(combine2["FC"].isna() | np.isinf(combine2["FC"]))
    ]

    # RMSE
    rmse1 = np.sqrt(((combine1["FC"] - combine1["log2FC"]) ** 2).mean())
    rmse2 = np.sqrt(((combine2["FC"] - combine2["log2FC"]) ** 2).mean())

    # Pearson (need >=2 points)
    r1 = scipy.stats.pearsonr(
        combine1["FC"], combine1["log2FC"]
    )[0] if len(combine1) > 1 else np.nan

    r2 = scipy.stats.pearsonr(
        combine2["FC"], combine2["log2FC"]
    )[0] if len(combine2) > 1 else np.nan

    return {
        "N_truth1": len(combine1),
        "N_truth2": len(combine2),
        "RMSE_truth1": rmse1,
        "RMSE_truth2": rmse2,
        "Pearson_truth1": r1,
        "Pearson_truth2": r2
    }


# --------------------------------------------------
# Main
# --------------------------------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Evaluate isoform FC accuracy against ground truth"
    )
    parser.add_argument("--expr", required=True, help="Expression matrix (tsv)")
    parser.add_argument("--truth", required=True, help="Quartet reference datasets")
    parser.add_argument("--truth2", required=True, help="RT-qPCR reference datasets")
    parser.add_argument("--lab", required=True, help="Lab prefix, e.g. lab1_")

    args = parser.parse_args()

    df = pd.read_csv(args.expr, sep="\t")
    truth = pd.read_csv(args.truth, sep="\t")
    truth2 = pd.read_csv(args.truth2, sep="\t")

    truth["transcript_id"] = truth["sample_pair"] + "_" + truth["transcript_id"]
    truth2["transcript_id"] = truth2["compare"] + "_" + truth2["transcript_id"]

    results = evaluate_fc_accuracy(
        df_subset=df,
        lab=args.lab,
        truth=truth,
        truth2=truth2
    )

    for k, v in results.items():
        print(f"{k}\t{v}")
