#!/usr/bin/env python3

import argparse
import logging
import os
import sys

# ===== Import internal modules from scripts/ =====
from scripts.unify import unify_results


def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(
        description="Unify splicing event results from multiple tools"
    )

    parser.add_argument(
        "--software", "-s",
        nargs="+",
        required=True,
        help="Tools used (e.g. SUPPA2 rMATS PSI-Sigma MAJIQ Spladder Whippet)"
    )

    parser.add_argument(
        "--event", "-e",
        nargs="+",
        required=False,
        choices=["SE", "A3SS", "A5SS", "MX", "RI", "AF", "AL"],
        help=(
            "Splicing event types:\n"
            "SE: Skipping Exon\n"
            "A3SS: Alternative 3' splice site\n"
            "A5SS: Alternative 5' splice site\n"
            "MX: Mutually Exclusive Exon\n"
            "RI: Retained Intron\n"
            "AF: Alternative First Exon\n"
            "AL: Alternative Last Exon"
        )
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input file or directory"
    )

    parser.add_argument(
        "--gtf", "-g",
        required=True,
        help="GTF annotation file"
    )

    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output file"
    )

    parser.add_argument(
        "--sample_name", "-sn",
        required=True,
        help="Sample name (must match folder name)"
    )

    parser.add_argument(
        "--groupA", "-A",
        nargs="+",
        required=True,
        help="space-seperated list of names of technical replicate for Sample group A, such as test_01 test_02 test_03"
    )

    parser.add_argument(
        "--groupB", "-B",
        nargs="+",
        required=True,
        help="space-seperated list of names of technical replicate for Sample group B, such as control_01 control_02 control_03"
    )

    return parser.parse_args()


def check_inputs(args):
    """
    Perform basic input validation
    """

    # Check input path
    if not os.path.exists(args.input):
        sys.exit(f"[ERROR] Input path not found: {args.input}")

    # Check GTF file
    if not os.path.exists(args.gtf):
        sys.exit(f"[ERROR] GTF file not found: {args.gtf}")

    # Check sample groups
    if len(args.groupA) == 0 or len(args.groupB) == 0:
        sys.exit("[ERROR] groupA and groupB must not be empty")


def main():
    """
    Main execution function
    """

    # Configure logging format and level
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    # Parse arguments
    args = parse_args()

    logging.info("Checking input parameters...")
    check_inputs(args)


    # Run main pipeline
    logging.info("Running unify_results...")
    unify_results(
        args.input,
        args.software,
        args.event,
        args.sample_name,
        args.output,
        args.gtf,
        args.groupA,
        args.groupB
    )

    logging.info("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
