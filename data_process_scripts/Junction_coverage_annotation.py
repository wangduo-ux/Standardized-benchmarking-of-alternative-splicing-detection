#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Transcript-level gene body coverage analysis (example script)

This script is modified from the geneBody_coverage.py script of RSeQC (v5.0.4)

This script demonstrates how to:
1. Compute gene body percentiles from a BED12 annotation
2. Calculate coverage along gene bodies using BAM files
3. Derive coverage uniformity metrics for a single sample group

This script is intended as an example for method illustration
and is suitable for sharing on GitHub.
"""

import os
import sys
import collections
import pysam
import pandas as pd
from time import strftime
from numpy import std, mean

# ============================================================
# Logging
# ============================================================
def printlog(msg):
    """Print progress message with timestamp"""
    msg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + msg
    print(msg, file=sys.stderr)


# ============================================================
# Gene body percentile calculation
# ============================================================
def genebody_percentile(refbed, min_mrna_len=100):
    """
    Calculate 100-percentile positions along gene bodies

    Parameters
    ----------
    refbed : str
        BED12 file describing transcript models
    min_mrna_len : int
        Minimum transcript length to keep

    Returns
    -------
    dict
        {transcript_id: (chrom, strand, percentile_positions)}
    """
    gene_percentiles = {}
    transcript_count = 0

    for line in open(refbed):
        if line.startswith(("#", "track", "browser")):
            continue

        try:
            fields = line.strip().split()
            chrom = fields[0]
            tx_start = int(fields[1])
            tx_end = int(fields[2])
            strand = fields[3]
            gene_name = fields[4]

            exon_sizes = list(map(int, fields[5].rstrip(",").split(",")))
            exon_starts = list(map(int, fields[6].rstrip(",").split(",")))
            exon_starts = [tx_start + x for x in exon_starts]
            exon_ends = [s + l for s, l in zip(exon_starts, exon_sizes)]

            transcript_count += 1
        except Exception:
            continue

        gene_bases = []
        for s, e in zip(exon_starts, exon_ends):
            gene_bases.extend(range(s + 1, e + 1))  # 1-based

        if len(gene_bases) < min_mrna_len:
            continue

        gene_id = "_".join(map(str, [chrom, tx_start, tx_end, gene_name, strand]))
        gene_percentiles[gene_id] = (chrom, strand, percentile_list(gene_bases))

    printlog(f"Loaded {transcript_count} transcripts")
    return gene_percentiles


def percentile_list(pos_list, n=100):
    """Generate n percentile positions from a list of coordinates"""
    return [pos_list[int(len(pos_list) * i / n)] for i in range(n)]


# ============================================================
# Gene body coverage
# ============================================================
def genebody_coverage(bam, gene_percentiles):
    """
    Calculate aggregated gene body coverage

    Parameters
    ----------
    bam : pysam.AlignmentFile
    gene_percentiles : dict

    Returns
    -------
    dict
        Aggregated coverage across percentiles
    """
    aggregated_cov = collections.defaultdict(int)

    for chrom, strand, positions in gene_percentiles.values():
        coverage = {p: 0 for p in positions}
        chrom_start = max(positions[0] - 1, 0)
        chrom_end = positions[-1]

        for col in bam.pileup(chrom, chrom_start, chrom_end, truncate=True):
            pos = col.pos + 1
            if pos not in coverage:
                continue

            count = 0
            for read in col.pileups:
                if read.is_del:
                    continue
                aln = read.alignment
                if aln.is_unmapped or aln.is_secondary or aln.is_duplicate:
                    continue
                count += 1

            coverage[pos] = count

        values = [coverage[p] for p in sorted(coverage)]
        if strand == "-":
            values = values[::-1]

        for i, v in enumerate(values):
            aggregated_cov[i] += v

    return aggregated_cov


# ============================================================
# Main example
# ============================================================
def main():

    # ---------------------------
    # User-defined inputs
    # ---------------------------
    BED_FILE = "/path/to/annotation.bed"
    BAM_FILE = "/path/to/example.bam"
    TRANSCRIPT_ID = "ENST00000XXXX"

    # ---------------------------
    # Load annotation
    # ---------------------------
    gene_percentiles = genebody_percentile(BED_FILE)

    # Keep only one transcript as an example
    gene_percentiles = {
        k: v for k, v in gene_percentiles.items()
        if TRANSCRIPT_ID in k
    }

    if not gene_percentiles:
        print("Transcript not found in annotation.")
        return

    # ---------------------------
    # Coverage calculation
    # ---------------------------
    bam = pysam.AlignmentFile(BAM_FILE, "rb")
    coverage = genebody_coverage(bam, gene_percentiles)

    df = pd.DataFrame(
        list(coverage.items()),
        columns=["percentile", "coverage"]
    )

    df["uniform"] = df["coverage"] / df["coverage"].max()
    uniform_prop = (df["uniform"] > 0.75).mean()

    print("Mean coverage:", df["coverage"].mean())
    print("Uniformity (>0.75):", uniform_prop)


if __name__ == "__main__":
    main()
