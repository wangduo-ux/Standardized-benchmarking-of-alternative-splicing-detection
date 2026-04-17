#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
STAR junction evaluation script (single lab, stdout output)

This script evaluates splice junction detection results from STAR
for ONE specified lab, using annotated and novel junction truth sets.

Results are printed directly to stdout instead of being saved to Excel.
"""

import argparse
import os
import pandas as pd


# =========================
# Argument parsing
# =========================
def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate STAR splice junction detection performance (single lab)"
    )
    parser.add_argument("--star_path", required=True, help="Path to STAR output directory")
    parser.add_argument("--annotation_file", required=True, help="SJDB annotation file, e.g. sjdbList.geneid.fromGTF.out.tab")
    parser.add_argument("--truth_annotated", required=True, help="Annotated junction truth file")
    parser.add_argument("--truth_novel", required=True, help="Novel junction truth file")
    parser.add_argument("--lab", required=True, help="Lab ID, e.g. lab17")
    return parser.parse_args()


# =========================
# Load junction annotation
# =========================
def load_annotation(annotation_file):
    """
    Load STAR sjdb annotation and generate junction IDs
    """
    annotate = pd.read_table(annotation_file, sep="\t", low_memory=False)
    annotate["strand"] = annotate["strand"].replace({"-": 2, "+": 1})
    annotate["ID_2"] = (
        annotate["seqnames"].astype(str)
        + "_"
        + annotate["start"].astype(str)
        + "_"
        + annotate["end"].astype(str)
        + "_"
        + annotate["strand"].astype(str)
    )
    annotate = annotate[["ID_2", "gene_id"]]
    annotate["anno"] = 1
    return annotate


# =========================
# Process one quartet
# =========================
def process(star_path, lab, samples, annotate):
    """
    Read STAR SJ.out.tab files for three replicates
    and merge junction read counts
    """
    dfs = []
    for s in samples:
        sj_file = os.path.join(star_path, f"{lab}_{s}.SJ.out.tab")
        df = pd.read_table(sj_file, header=None, sep="\t", low_memory=False)
        df.columns = [
            "seqnames", "start", "end", "strand",
            "motif", "annotate", "unique_reads",
            "multi_reads", "splicing_overhang"
        ]
        df["reads"] = df["unique_reads"] + df["multi_reads"]
        df["ID_2"] = (
            df["seqnames"].astype(str)
            + "_"
            + df["start"].astype(str)
            + "_"
            + df["end"].astype(str)
            + "_"
            + df["strand"].astype(str)
        )
        dfs.append(df[["ID_2", "reads"]])

    merged = pd.concat(dfs, axis=1)
    merged.columns = ["jun1", "jun2", "jun3"]

    # Keep junctions with at least one supporting read
    merged = merged.loc[merged.sum(axis=1) > 0].reset_index()

    merged = pd.merge(merged, annotate, on="ID_2", how="left")
    merged["anno"] = merged["anno"].fillna(0).astype(int)

    annotated = merged[merged["anno"] == 1].copy()
    novel = merged[merged["anno"] == 0].copy()

    # Remove strand information for downstream matching
    annotated["ID_2"] = annotated["ID_2"].str.split("_").str[:3].str.join("_")
    novel["ID_2"] = novel["ID_2"].str.split("_").str[:3].str.join("_")

    return annotated, novel


# =========================
# Combine quartet results
# =========================
def concat(star_path, lab, samples, quartet, annotate, truth):
    """
    Add quartet prefix and filter by gene-level truth
    """
    annotated, novel = process(star_path, lab, samples, annotate)

    annotated["status_test"] = (
        annotated[["jun1", "jun2", "jun3"]]
        .ge(1).sum(axis=1).ge(2)
        .map({True: "Y_H", False: "N_H"})
    )
    annotated["ID_2"] = quartet + "_" + annotated["ID_2"]

    truth_q = truth[truth["ID_2"].str.startswith(quartet)]
    annotated = annotated[annotated["gene_id"].isin(truth_q["gene_id"])]

    novel["status"] = (
        novel[["jun1", "jun2", "jun3"]]
        .ge(1).sum(axis=1).ge(2)
        .map({True: "Y_H", False: "N_H"})
    )
    novel = novel[novel["status"] == "Y_H"]
    novel["ID_2"] = quartet + "_" + novel["ID_2"]

    return annotated, novel


# =========================
# Classification logic
# =========================
def classify(row):
    """
    Assign TP / FP / FN / TN
    """
    if row["status"] == "Y_H" and row["status_test"] == "Y_H":
        return "TP"
    elif row["status"] == "N_H" and row["status_test"] == "N_H":
        return "TN"
    elif pd.isna(row["status"]) and row["status_test"] == "N_H":
        return "TN"
    elif pd.isna(row["status"]) and row["status_test"] == "Y_H":
        return "FP"
    elif row["status"] == "Y_H" and row["status_test"] == "N_H":
        return "FN"
    elif row["status"] == "N_H" and row["status_test"] == "Y_H":
        return "FP"


# =========================
# Main workflow
# =========================
def main():
    args = parse_args()

    annotate = load_annotation(args.annotation_file)

    # Load annotated truth
    truth = pd.read_table(args.truth_annotated, sep="\t", low_memory=False)
    truth["ID_2"] = (
        truth["sample"] + "_"
        + truth["seqnames"].astype(str)
        + "_"
        + truth["start"].astype(str)
        + "_"
        + truth["end"].astype(str)
    )
    truth = truth[["ID_2", "status", "gene_id"]]

    # Load novel truth
    truth_novel = pd.read_table(args.truth_novel, sep="\t", low_memory=False)
    truth_novel["ID_2"] = truth_novel["sample"] + "_" + truth_novel["ID"]

    annotated_all = []
    novel_all = []

    for quartet, samples in zip(
        ["M8", "F7", "D5", "D6"],
        [
            [202207, 202208, 202209],
            [202210, 202211, 202212],
            [202213, 202214, 202215],
            [202222, 202223, 202224],
        ],
    ):
        a, n = concat(
            args.star_path,
            args.lab,
            samples,
            quartet,
            annotate,
            truth,
        )
        annotated_all.append(a)
        novel_all.append(n)

    annotated = pd.concat(annotated_all, ignore_index=True)
    novel = pd.concat(novel_all, ignore_index=True)
    novel["classification"] = "TP"

    novel_recall = (
        len(set(novel["ID_2"]) & set(truth_novel["ID_2"]))
        / len(truth_novel)
        if len(truth_novel) > 0 else 0
    )

    merged = pd.merge(annotated, truth, on="ID_2", how="outer")
    merged["status_test"] = merged["status_test"].fillna("N_H")
    merged["classification"] = merged.apply(classify, axis=1)

    TP = (merged["classification"] == "TP").sum()
    FP = (merged["classification"] == "FP").sum()
    FN = (merged["classification"] == "FN").sum()
    TN = (merged["classification"] == "TN").sum()

    # ===== stdout output =====
    print(
        "\t".join(
            map(
                str,
                [
                    args.lab,
                    TP,
                    FP,
                    FN,
                    TN,
                    round(novel_recall, 4),
                ],
            )
        )
    )


if __name__ == "__main__":
    main()
