#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import scipy.stats

def process_func(df,sample):
  df = df.copy()
  df.iloc[:, 1:7]= np.log2(df.iloc[:, 1:7] + 0.01)
  df = df.reset_index(drop=True)
  df['mean'] = df.iloc[:, 4:7].mean(axis=1)
  df.iloc[:, 1:4] = df.iloc[:, 1:4].sub(df.iloc[:, -1], axis=0)
  df['mean'] = df.iloc[:, 1:4].mean(axis=1)
  data = df.loc[:, ["transcript_id", "mean"]]
  data = data.rename(columns={'mean': sample})
  return data


def calculate_mix_accuracy(df, lab):
    """
    Calculate correlation and RMSE between observed and expected expression
    for T1 / T2 samples based on M8 reference.

    Parameters
    ----------
    df : pandas.DataFrame
        Input expression matrix
    selected_cols : list
        Columns to select from df
    lab : str
        Lab prefix, e.g. "lab1_"
    process_func : function
        Function to process expression matrices (user-defined)

    Returns
    -------
    dict
        Pearson correlation and RMSE for T1 and T2
    """

    # -----------------------------
    # Calculate mean expression (optional, kept from original code)
    # -----------------------------
    test = df.copy()

    # -----------------------------
    # Prepare datasets
    # -----------------------------
    data1 = df[
        ["transcript_id",
         f"{lab}202207", f"{lab}202208", f"{lab}202209",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    data4 = df[
        ["transcript_id",
         f"{lab}202216", f"{lab}202217", f"{lab}202218",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    data5 = df[
        ["transcript_id",
         f"{lab}202219", f"{lab}202220", f"{lab}202221",
         f"{lab}202222", f"{lab}202223", f"{lab}202224"]
    ]

    # -----------------------------
    # Expression filtering
    # -----------------------------
    data1 = data1[data1.iloc[:, 1:7].mean(axis=1) > 0]
    data4 = data4[data4.iloc[:, 1:7].mean(axis=1) > 0]
    data5 = data5[data5.iloc[:, 1:7].mean(axis=1) > 0]

    # -----------------------------
    # Process data
    # -----------------------------
    data1 = process_func(data1, "M8")
    data4 = process_func(data4, "T1")
    data5 = process_func(data5, "T2")

    # -----------------------------
    # Merge datasets
    # -----------------------------
    merged = pd.merge(data1, data4, on="transcript_id")
    merged = pd.merge(merged, data5, on="transcript_id")

    dff = merged[["transcript_id", "M8", "T1", "T2"]].reset_index(drop=True)

    # -----------------------------
    # Expected expression (mixture model)
    # -----------------------------
    dff["expected_T1"] = np.log2(
        0.245182377 + 0.745119679 * (2 ** dff["M8"])
    )

    dff["expected_T2"] = np.log2(
        0.739998388 + 0.240258211 * (2 ** dff["M8"])
    )

    # -----------------------------
    # Metrics
    # -----------------------------
    corr_T1, _ = scipy.stats.pearsonr(dff["T1"], dff["expected_T1"])
    corr_T2, _ = scipy.stats.pearsonr(dff["T2"], dff["expected_T2"])

    rmse_T1 = np.sqrt(((dff["T1"] - dff["expected_T1"]) ** 2).mean())
    rmse_T2 = np.sqrt(((dff["T2"] - dff["expected_T2"]) ** 2).mean())

    return {
        "Pearson_T1": corr_T1,
        "Pearson_T2": corr_T2,
        "RMSE_T1": rmse_T1,
        "RMSE_T2": rmse_T2,
        "N_transcripts": dff.shape[0]
    }


# ============================================================
# Example usage
# ============================================================
if __name__ == "__main__":

    import argparse
    import pandas as pd

    parser = argparse.ArgumentParser(
        description="Evaluate expression accuracy based on the revovery of mixing ratios"
    )
    parser.add_argument(
        "--expr",
        required=True,
        help="The file path for the data is represented, with columns for transcript_id and sample number starting with lab."
    )
    parser.add_argument(
        "--lab",
        required=True,
        help="Lab prefix, e.g. lab1_"
    )

    args = parser.parse_args()

    # -----------------------------
    # Read expression matrix
    # -----------------------------
    df = pd.read_csv(
        args.expr,
        sep="\t",
        header=0,
        dtype={0: str}
    )

    if "transcript_id" not in df.columns:
        raise ValueError("Input expression matrix must contain 'transcript_id' column")


    # -----------------------------
    # Run evaluation
    # -----------------------------
    results = calculate_mix_accuracy(
        df=df,
        selected_cols=selected_cols,
        lab=args.lab,
        process_func=pr_
