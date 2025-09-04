# scripts/__init__.py

# Expose main functions from each module for easy import
from .file_loader import load_file
from .file_io import save_dataframe


from scripts.registry import (
    PARSERS,
    get_support_software_list,
    parse_suppa2,
    parse_rmats,
    parse_majiq,
    parse_psi_sigma
)

# --- File I/O ---
from .file_loader import load_file
from .file_io import save_dataframe

# --- DAS / ASE processing ---
from .das_parser import annotate_das_class,ase_extract
from .das_plot_utils import (
    plot_venn_per_event_vennlib,
    plot_upset_per_event_combined,
    plot_multiple_event_summaries
)

# --- dPSI plotting & extraction ---
from .dpsi_plot_utils import plot_all_events_grid, extract_dpsi, plot_corr_dpsi

# --- rMATS novel filtering ---
from .filter_novel_rmats import filter_rmats_novel_events

# --- GTF parsing ---
from .gtf_parser import extract_gene_transcript_strand

from scripts.suppa2_parser import parse_SUPPA2
from scripts.rmats_parser import parse_rMATS_SE, parse_rMATS_A3SS_A5SS, parse_rMATS_RI, parse_rMATS_MX
from scripts.majiq_parser import (
    parse_MAJIQ_SE, parse_MAJIQ_A3SS_A5SS, parse_MAJIQ_AF_AL,
    parse_MAJIQ_RI, parse_MAJIQ_MX
)
from scripts.psisigma_parser import (
    parse_psisigma_SE_A3SS_A5SS, parse_psisigma_RI, parse_psisigma_MX
)

