import re
import pandas as pd

def parse_rMATS_SE(row, event_type):
    """
    Parse rMATS Skipped Exon (SE) results to construct a uniform_ID.

    Parameters
    ----------
    row : pd.Series or list
        A single row of rMATS SE results.
    event_type : str
        The event type, expected to be "SE".

    Returns
    -------
    str or None
        A uniform event ID string, or None if parsing fails.
    """
    (ID, GeneID, geneSymbol, Chr, strand,
     exonStart_0base, exonEnd, upstreamES, upstreamEE,
     downstreamES, downstreamEE,
     ID2, IC_SAMPLE_1, SC_SAMPLE_1, IC_SAMPLE_2, SC_SAMPLE_2,
     IncFormLen, SkipFormLen, PValue, FDR,
     IncLevel1, IncLevel2, IncLevelDifference) = row

    m = re.match(r'^"(ENSG[^"]+)"$', str(row.GeneID))
    gene_id = m.group(1) if m else str(row.GeneID).strip('"')

    Chr = str(Chr)
    if Chr.lower().startswith("chr"):
        Chr = Chr[3:]

    try:
        up_ee = int(row.upstreamEE)
        exon0 = int(row.exonStart_0base)
        exon_e = int(row.exonEnd)
        down_es = int(row.downstreamES)
    except ValueError:
        return None

    start1 = exon0 + 1
    down1 = down_es + 1
    return f"{gene_id};{event_type}:{Chr}:{up_ee}-{start1}:{exon_e}-{down1}:{row.strand}"


def parse_rMATS_A3SS_A5SS(row, event_type):
    """
    Parse rMATS A3SS or A5SS results to construct a uniform_ID.

    Parameters
    ----------
    row : pd.Series or list
        A single row of rMATS A3SS or A5SS results.
    event_type : str
        The event type, "A3SS" or "A5SS".

    Returns
    -------
    str or None
        A uniform event ID string, or None if parsing fails.
    """
    (ID, GeneID, geneSymbol, Chr, strand,
     longExonStart_0base, longExonEnd, shortES, shortEE,
     flankingES, flankingEE,
     ID2, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2,
     IncFormLen, SkipFormLen, PValue, FDR,
     IncLevel1, IncLevel2, IncLevelDifference) = row

    m = re.match(r'^"(ENSG[^"]+)"$', str(GeneID))
    gene_id = m.group(1) if m else str(GeneID).strip('"')

    Chr = str(Chr)
    if Chr.lower().startswith("chr"):
        Chr = Chr[3:]

    try:
        longExonStart_0base = int(longExonStart_0base)
        longExonEnd = int(longExonEnd)
        shortES = int(shortES)
        shortEE = int(shortEE)
        flankingES = int(flankingES)
        flankingEE = int(flankingEE)
    except ValueError:
        return None

    if event_type == "A3SS":
        if strand == "+":
            return f"{gene_id};A3:{Chr}:{flankingEE}-{longExonStart_0base+1}:{flankingEE}-{shortES+1}:{strand}"
        else:
            return f"{gene_id};A3:{Chr}:{longExonEnd}-{flankingES+1}:{shortEE}-{flankingES+1}:{strand}"

    elif event_type == "A5SS":
        if strand == "+":
            return f"{gene_id};A5:{Chr}:{longExonEnd}-{flankingES+1}:{shortEE}-{flankingES+1}:{strand}"
        else:
            return f"{gene_id};A5:{Chr}:{flankingEE}-{longExonStart_0base+1}:{flankingEE}-{shortES+1}:{strand}"

    return None


def parse_rMATS_RI(row, event_type) :
    """
    Parse rMATS Retained Intron (RI) results to construct a uniform_ID.

    Parameters
    ----------
    row : pd.Series or list
        A single row of rMATS RI results.
    event_type : str
        The event type, expected to be "RI".

    Returns
    -------
    str or None
        A uniform event ID string, or None if parsing fails.
    """
    (ID, GeneID, geneSymbol, Chr, strand,
     riExonStart_0base, riExonEnd, upstreamES, upstreamEE,
     downstreamES, downstreamEE,
     ID2, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2,
     IncFormLen, SkipFormLen, PValue, FDR,
     IncLevel1, IncLevel2, IncLevelDifference) = row

    m = re.match(r'^"(ENSG[^"]+)"$', str(GeneID))
    gene_id = m.group(1) if m else str(GeneID).strip('"')

    Chr = str(Chr)
    if Chr.lower().startswith("chr"):
        Chr = Chr[3:]

    try:
        upstreamES = int(upstreamES)
        upstreamEE = int(upstreamEE)
        downstreamES = int(downstreamES)
        downstreamEE = int(downstreamEE)
    except ValueError:
        return None

    return f"{gene_id};RI:{Chr}:{upstreamES+1}:{upstreamEE}-{downstreamES+1}:{downstreamEE}:{strand}"


def parse_rMATS_MX(row, event_type) :
    """
    Parse rMATS Mutually Exclusive Exon (MX) results to construct a uniform_ID.

    Parameters
    ----------
    row : pd.Series or list
        A single row of rMATS MX results.
    event_type : str
        The event type, expected to be "MX".

    Returns
    -------
    str or None
        A uniform event ID string, or None if parsing fails.
    """
    (ID, GeneID, geneSymbol, Chr, strand,
     exon1_start, exon1_end, exon2_start, exon2_end,
     upstreamES, upstreamEE, downstreamES, downstreamEE,
     ID2, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2,
     IncFormLen, SkipFormLen, PValue, FDR,
     IncLevel1, IncLevel2, IncLevelDifference) = row

    m = re.match(r'^"(ENSG[^"]+)"$', str(GeneID))
    gene_id = m.group(1) if m else str(GeneID).strip('"')

    Chr = str(Chr)
    if Chr.lower().startswith("chr"):
        Chr = Chr[3:]

    try:
        exon1_start = int(exon1_start)
        exon1_end = int(exon1_end)
        exon2_start = int(exon2_start)
        exon2_end = int(exon2_end)
        upstreamEE = int(upstreamEE)
        downstreamES = int(downstreamES)
    except ValueError:
        return None

    return (
        f"{gene_id};MX:{Chr}:{upstreamEE}-{exon1_start+1}:{exon1_end}-{downstreamES+1}:"
        f"{upstreamEE}-{exon2_start+1}:{exon2_end}-{downstreamES+1}:{strand}"
    )
