"""
ASI.py

Module functionality:
- Load alternative splicing (AS) event data from a specified directory
- Supports software such as SUPPA2, MAJIQ, rMATS, PSI-Sigma
- Extract DSE for each sample and event type
- Returns intergration results for downstream analyses

Dependencies:
- pandas
- os
- Custom functions: annotate_das_class
"""

import os
import pandas as pd
import argparse
from scripts import (
    load_files,
    integration_analysis,
    annotate_dse,
    select_rows
)
  

def main():
    parser = argparse.ArgumentParser(
        description="Unify splicing event results from multiple tools"
    )
    parser.add_argument('--software', '-s',
                        nargs='+',
                        #choices=PARSERS.keys(),
                        required=True,
                        help='space separated list of tools for benchmarking from the following list:SUPPA2, rMATS, PSI-Sigma, MAJIQ')
    parser.add_argument("-e", "--event", nargs='+', required=False, choices=["SE", "A3SS", "A5SS","MX", "RI", "AF", "AL"],
                        help="list of events to analyze. "
                        "(space separated)\n\n"
                        "Options:\n"
                        "\tSE -- Skipping Exon\n"
                        "\tA3SS -- Alternative Splice Site (3')\n"
                        "\tA5SS -- Alternative Splice Site (5')\n"
                        "\tMX -- Mutually Exclusive Exon\n"
                        "\tRI -- Retained Intron\n"
                        "\tAF -- Alternative First Exon\n"
                        "\tAL -- Alternative Last Exon\n")
    parser.add_argument('--input', '-i',
                        required=True,
                        help='Input file or directory')
    parser.add_argument('--output', '-o',
                        help='Output file')
    parser.add_argument('--sample_name', '-sn',
                        required=True,
                        help='The sample name needs to correspond to the folder name')
    args = parser.parse_args()
    integration_analysis(args.input, args.software, args.event, args.sample_name, args.output)


if __name__ == '__main__':
    main()
