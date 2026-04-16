"""
Registry of parser functions for different splicing software.
"""
# ===== imports =====

from scripts.suppa2_parser import parse_SUPPA2
from scripts.rmats_parser import (
    parse_rMATS_SE,
    parse_rMATS_A3SS_A5SS,
    parse_rMATS_RI,
    parse_rMATS_MX
)
from scripts.majiq_parser import (
    parse_MAJIQ_SE,
    parse_MAJIQ_A3SS_A5SS,
    parse_MAJIQ_AF_AL,
    parse_MAJIQ_RI,
    parse_MAJIQ_MX
)
from scripts.psisigma_parser import (
    parse_psisigma_SE_A3SS_A5SS,
    parse_psisigma_RI,
    parse_psisigma_MX
)
from scripts.spladder_parser import (
    parse_Spladder_SE,
    parse_Spladder_A3SS_A5SS,
    parse_Spladder_RI,
    parse_Spladder_MX
)
from scripts.whippet_parser import parse_Whippet


import re

# ===== software name normalization =====
def normalize_software_name(name):
    name = name.lower()

    if re.search("rmats", name):
        return "rMATS"

    if re.search("majiq", name):
        return "MAJIQ"

    if re.search("spladder", name):
        return "Spladder"

    if re.search("suppa", name):
        return "SUPPA2"

    if re.search("whippet", name):
        return "Whippet"

    if re.search("sigma", name):
        return "PSI-Sigma"

    return name


# ===== parser mapping =====
PARSER_FUNCTIONS = {
    'SUPPA2': {
        'SE': parse_SUPPA2,
        'A5SS': parse_SUPPA2,
        'A3SS': parse_SUPPA2,
        'AF': parse_SUPPA2,
        'AL': parse_SUPPA2,
        'RI': parse_SUPPA2,
        'MX': parse_SUPPA2
    },
    'rMATS': {
        'SE': parse_rMATS_SE,
        'A5SS': parse_rMATS_A3SS_A5SS,
        'A3SS': parse_rMATS_A3SS_A5SS,
        'RI': parse_rMATS_RI,
        'MX': parse_rMATS_MX
    },
    'MAJIQ': {
        'SE': parse_MAJIQ_SE,
        'A5SS': parse_MAJIQ_A3SS_A5SS,
        'A3SS': parse_MAJIQ_A3SS_A5SS,
        'AF': parse_MAJIQ_AF_AL,
        'AL': parse_MAJIQ_AF_AL,
        'RI': parse_MAJIQ_RI,
        'MX': parse_MAJIQ_MX
    },
    'PSI-Sigma': {
        'SE': parse_psisigma_SE_A3SS_A5SS,
        'A5SS': parse_psisigma_SE_A3SS_A5SS,
        'A3SS': parse_psisigma_SE_A3SS_A5SS,
        'RI': parse_psisigma_RI,
        'MX': parse_psisigma_MX
    },
    'Spladder': {
        'SE': parse_Spladder_SE,
        'A5SS': parse_Spladder_A3SS_A5SS,
        'A3SS': parse_Spladder_A3SS_A5SS,
        'RI': parse_Spladder_RI,
        'MX': parse_Spladder_MX
    },
    'Whippet': {
        'SE': parse_Whippet,
        'A5SS': parse_Whippet,
        'A3SS': parse_Whippet,
        'RI': parse_Whippet,
        'AF': parse_Whippet,
        'AL': parse_Whippet
    }
}


# ===== main dispatcher =====
def parse_software(df, software, event_type):
    """
    Dispatch parsing function based on software and event type.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame
    software : str
        Software name
    event_type : str
        Splicing event type

    Returns
    -------
    pd.DataFrame
        Parsed and standardized DataFrame
    """

    software = normalize_software_name(software)

    if software not in PARSER_FUNCTIONS:
        raise ValueError(f"Unsupported software: {software}")

    if event_type not in PARSER_FUNCTIONS[software]:
        raise ValueError(
            f"Unsupported event type '{event_type}' for {software}"
        )

    parser = PARSER_FUNCTIONS[software][event_type]

    if parser is None:
        raise ValueError(f"No parser defined for {software} {event_type}")

    return parser(df, event_type)
