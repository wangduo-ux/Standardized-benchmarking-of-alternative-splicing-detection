# =========================
# package init
# =========================

# IO
from .io import load_file

# ASE
from .ase_parser import  ase_extract

from .dpsi_parser import extract_dpsi
# parser
from .parser import (
    parse_SUPPA2,

    parse_rMATS_SE,
    parse_rMATS_A3SS_A5SS,
    parse_rMATS_RI,
    parse_rMATS_MX,

    parse_psisigma_SE_A3SS_A5SS,
    parse_psisigma_MX,
    parse_psisigma_RI,

    parse_MAJIQ_SE,
    parse_MAJIQ_A3SS_A5SS,
    parse_MAJIQ_AF_AL,
    parse_MAJIQ_RI,
    parse_MAJIQ_MX,

    parse_Spladder_SE,
    parse_Spladder_A3SS_A5SS,
    parse_Spladder_RI,
    parse_Spladder_MX,

    parse_Whippet
)

# pairwise / merge
from .pairwise_analysis import (
    build_pairwise_df,
    compute_pairwise_pearson,
    summarize_whippet_overlap
)

from .event_merger import merge_splicing_events
    
    
__all__ = [
    "load_file",

    "extract_dpsi",
    "ase_extract",

    "parse_SUPPA2",

    "parse_rMATS_SE",
    "parse_rMATS_A3SS_A5SS",
    "parse_rMATS_RI",
    "parse_rMATS_MX",

    "parse_psisigma_SE_A3SS_A5SS",
    "parse_psisigma_MX",
    "parse_psisigma_RI",

    "parse_MAJIQ_SE",
    "parse_MAJIQ_A3SS_A5SS",
    "parse_MAJIQ_AF_AL",
    "parse_MAJIQ_RI",
    "parse_MAJIQ_MX",

    "parse_Spladder_SE",
    "parse_Spladder_A3SS_A5SS",
    "parse_Spladder_RI",
    "parse_Spladder_MX",

    "parse_Whippet",

    "assign_event_id",
    "fix_majiq_na",
    "merge_splicing_events",
    "build_pairwise_df",
    "compute_pairwise_pearson",
    "normalize_software_name",
    "build_upset_df",
    "summarize_whippet_overlap"
]
