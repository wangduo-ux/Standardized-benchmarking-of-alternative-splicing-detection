import re
import pandas as pd
import itertools

def parse_psisigma_SE_A3SS_A5SS(df, event_type):
    """
    Parse SE, A3SS, and A5SS results output by PSI-Sigma and construct uniform_IDs.
    """

    def generate_ID(row):
        (
            Event_Region, Gene_Symbol, Target_Exon, Event_Type, N, T, Exon_Type,
            Reference_Transcript, ΔPSI, p_value, FDR, N_Values, T_Values,
            Database_ID, type, gene_id, strand
        ) = row

        # Normalize gene_id
        m = re.match(r'^"(ENSG[^"]+)"$', str(row.gene_id))
        gene_id = m.group(1) if m else str(row.gene_id).strip('"')

        try:
            # Parse target exon coordinates
            Target_Exon_coord = re.match(r'^([^:]+):(\d+)-(\d+)$', Target_Exon)
            Chr = str(Target_Exon_coord.group(1))
            if Chr.lower().startswith('chr'):
                Chr = Chr[3:]
            exon_start = int(Target_Exon_coord.group(2))
            exon_end = int(Target_Exon_coord.group(3))

            # Parse event region coordinates
            Event_Region_coord = re.match(r'^([^:]+):(\d+)-(\d+)$', Event_Region)
            up_ee = int(Event_Region_coord.group(2))
            down_es = int(Event_Region_coord.group(3))
        except ValueError:
            return None

        # Construct uniform_ID depending on event type
        if Event_Type == "SES":
            up_ee = up_ee - 1
            down_es = down_es + 1
            uniform_id = (
                f"{gene_id};SE:{Chr}:{up_ee}-{exon_start}:"
                f"{exon_end}-{down_es}:{row.strand}"
            )
        elif Event_Type in ["A3SS", "TSS|A3SS"]:
            if strand == "+":
                uniform_id = (
                    f"{gene_id};A3:{Chr}:{up_ee-1}-{exon_start}:"
                    f"{up_ee-1}-{down_es+1}:{strand}"
                )
            elif strand == "-":
                uniform_id = (
                    f"{gene_id};A3:{Chr}:{exon_end}-{down_es+1}:"
                    f"{up_ee-1}-{down_es+1}:{strand}"
                )
        elif Event_Type in ["A5SS", "TSS|A5SS"]:
            if strand == "+":
                uniform_id = (
                    f"{gene_id};A5:{Chr}:{exon_end}-{down_es+1}:"
                    f"{up_ee-1}-{down_es+1}:{strand}"
                )
            elif strand == "-":
                uniform_id = (
                    f"{gene_id};A5:{Chr}:{up_ee-1}-{exon_start}:"
                    f"{up_ee-1}-{down_es+1}:{strand}"
                )
        return uniform_id

    # Filter by event type
    if event_type == "SE":
        df = df[df['Event Type'] == 'SES'].copy()
    elif event_type == "A3SS":
        df = df[df['Event Type'].isin(['A3SS', 'TSS|A3SS'])].copy()
    elif event_type == "A5SS":
        df = df[df['Event Type'].isin(['A5SS', 'TSS|A5SS'])].copy()

    # Apply uniform_ID generation
    df['uniform_ID'] = df.apply(generate_ID, axis=1)
    return df


def parse_psisigma_MX(df, event_type):
    """
    Parse MX (mutually exclusive exons) results output by PSI-Sigma and construct uniform_IDs.
    """
    df = df[df["Event Type"] == "MXS"]
    df = df[df['Event Region'].map(df['Event Region'].value_counts()) > 1]

    # Pre-allocate uniform_ID column
    df['uniform_ID'] = pd.NA
    merged_dfs = []

    for event_region, group in df.groupby("Event Region"):
        chrom, coord = event_region.split(":")
        region_start, region_end = map(int, coord.split('-'))
        gene_id = group["gene_id"].iloc[0]
        strand = group["strand"].iloc[0]

        exon_info = []
        for idx, exon_str in zip(group.index, group["Target Exon"]):
            try:
                _, coord_part = exon_str.split(":")
                start, end = map(int, coord_part.split('-'))
                exon_info.append((idx, start, end))
            except:
                continue

        exon_info = sorted(exon_info, key=lambda x: x[1])

        # Compare all exon pairs
        for (i1, s1, e1), (i2, s2, e2) in itertools.combinations(exon_info, 2):
            if e1 < s2 or e2 < s1:  # non-overlapping
                if s1 < s2:
                    left = f"{s1}:{e1}"
                    right = f"{s2}:{e2}"
                    keep_index = i1
                else:
                    left = f"{s2}:{e2}"
                    right = f"{s1}:{e1}"
                    keep_index = i2

                uniform_id = (
                    f"{gene_id};MX:{chrom}:{region_start-1}-{left}-{region_end+1}:"
                    f"{region_start-1}-{right}-{region_end+1}:{strand}"
                )

                merged_row = pd.concat([df.loc[i1], df.loc[i2]], axis=1).T
                merged_row["uniform_ID"] = uniform_id
                merged_row.reset_index(drop=True)
                merged_dfs.append(merged_row)

    df = pd.concat(merged_dfs, ignore_index=True)
    return df


def parse_psisigma_RI(df, event_type):
    """
    Parse RI (retained intron) results output by PSI-Sigma and construct uniform_IDs.
    """
    df = df[df['Event Type'].isin(['RI', 'IR (overlapping region)'])].copy()

    def generate_ID(row):
        (
            Event_Region, Gene_Symbol, Target_Exon, Event_Type, N, T, Exon_Type,
            Reference_Transcript, ΔPSI, p_value, FDR, N_Values, T_Values,
            Database_ID, type, gene_id, strand
        ) = row

        # Normalize gene_id
        m = re.match(r'^"(ENSG[^"]+)"$', str(row.gene_id))
        gene_id = m.group(1) if m else str(row.gene_id).strip('"')

        try:
            # Parse exon coordinates
            Target_Exon_coord = re.match(r'^([^:]+):(\d+)-(\d+)$', Target_Exon)
            Chr = str(Target_Exon_coord.group(1))
            if Chr.lower().startswith('chr'):
                Chr = Chr[3:]
            exon_start = int(Target_Exon_coord.group(2))
            exon_end = int(Target_Exon_coord.group(3))

            # Parse region coordinates
            Event_Region_coord = re.match(r'^([^:]+):(\d+)-(\d+)$', Event_Region)
            up_ee = int(Event_Region_coord.group(2))
            down_es = int(Event_Region_coord.group(3))
        except ValueError:
            return None

        # Construct uniform_ID
        uniform_id = (
            f"{gene_id};RI:{Chr}:{up_ee}:{exon_start-1}-{exon_end+1}:{down_es}:{row.strand}"
        )
        return uniform_id

    # Apply uniform_ID generation
    df['uniform_ID'] = df.apply(generate_ID, axis=1)
    return df
