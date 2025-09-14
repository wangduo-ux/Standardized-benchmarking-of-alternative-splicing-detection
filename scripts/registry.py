"""
Registry of parser functions for different splicing software.
"""

from scripts.suppa2_parser import parse_SUPPA2
from scripts.rmats_parser import parse_rMATS_SE, parse_rMATS_A3SS_A5SS, parse_rMATS_RI, parse_rMATS_MX
from scripts.majiq_parser import (
    parse_MAJIQ_SE, parse_MAJIQ_A3SS_A5SS, parse_MAJIQ_AF_AL,
    parse_MAJIQ_RI, parse_MAJIQ_MX
)
from scripts.psisigma_parser import (
    parse_psisigma_SE_A3SS_A5SS, parse_psisigma_RI, parse_psisigma_MX
)
from scripts.gtf_parser import extract_gene_transcript_strand
import pandas as pd


# Low-level mapping: software → event_type → parser
PARSER_FUNCTIONS = {
    'SUPPA2': {
        'SE': parse_SUPPA2, 'A5SS': parse_SUPPA2, 'A3SS': parse_SUPPA2,
        'AF': parse_SUPPA2, 'AL': parse_SUPPA2, 'RI': parse_SUPPA2, 'MX': parse_SUPPA2
    },
    'rMATS': {
        'SE': parse_rMATS_SE, 'A5SS': parse_rMATS_A3SS_A5SS, 'A3SS': parse_rMATS_A3SS_A5SS,
        'RI': parse_rMATS_RI, 'MX': parse_rMATS_MX
    },
    'MAJIQ': {
        'SE': parse_MAJIQ_SE, 'A5SS': parse_MAJIQ_A3SS_A5SS, 'A3SS': parse_MAJIQ_A3SS_A5SS,
        'AF': parse_MAJIQ_AF_AL, 'AL': parse_MAJIQ_AF_AL,
        'RI': parse_MAJIQ_RI, 'MX': parse_MAJIQ_MX
    },
    'PSI-Sigma': {
        'SE': parse_psisigma_SE_A3SS_A5SS, 'A5SS': parse_psisigma_SE_A3SS_A5SS,
        'A3SS': parse_psisigma_SE_A3SS_A5SS,
        'RI': parse_psisigma_RI, 'MX': parse_psisigma_MX
    }
}


# High-level parsing wrappers
def parse_suppa2(df, event_type):
    parser = PARSER_FUNCTIONS["SUPPA2"][event_type]
    return parser(df, event_type)


def parse_rmats(df, event_type):
    parser = PARSER_FUNCTIONS["rMATS"][event_type]
    df['uniform_ID'] = df.apply(lambda row: parser(row, event_type), axis=1)
    return df


def parse_majiq(df, event_type):
    parser = PARSER_FUNCTIONS["MAJIQ"][event_type]
    return parser(df, event_type)


def parse_psi_sigma(df, event_type, gtf_path):
    df = df.rename(columns={'Reference Transcript': 'transcript_id'})
    df['type'] = df['transcript_id'].astype(str).str.split('.').str[0].apply(lambda x: 'novel' if x == 'Ex' else '-')
    df['transcript_id'] = df['transcript_id'].astype(str).str.split('.').str[-1]
    gtf_df = extract_gene_transcript_strand(gtf_path)
    df = df.merge(gtf_df, on="transcript_id", how="left")
    parser = PARSER_FUNCTIONS["PSI-Sigma"][event_type]
    return parser(df, event_type)


# Master dictionary for external calls
PARSERS = {
    'SUPPA2': parse_suppa2,
    'rMATS': parse_rmats,
    'MAJIQ': parse_majiq,
    'PSI-Sigma': parse_psi_sigma,
}


def get_support_software_list(event_type):
    """Return list of software supporting a given event type."""
    if event_type in ('SE', 'A3SS', 'A5SS', 'MX', 'RI'):
        return ['SUPPA2', 'rMATS', 'PSI-Sigma', 'MAJIQ']
    elif event_type in ('AF', 'AL'):
        return ['SUPPA2', 'MAJIQ']
    else:
        return []
