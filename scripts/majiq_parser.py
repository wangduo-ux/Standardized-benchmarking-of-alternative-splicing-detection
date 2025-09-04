import re
import pandas as pd


def parse_MAJIQ_SE(df, event_type):
    """
    Parse MAJIQ SE (Skipped Exon) event DataFrame with 4 lines per event_id
    to build uniform_ID.
    """
    records = []
    for event_id, group in df.groupby('event_id'):
        group = group.reset_index(drop=True)
        raw_gene_id = str(group.loc[0, 'gene_id'])
        match = re.search(r'ENSG\d+', raw_gene_id)
        gene_id = match.group(0) if match else raw_gene_id
        strand = group.loc[0, 'strand']
        chrom = group.loc[0, 'seqid'].replace('seqid', '')

        try:
            c1c2_row = group[group['junction_name'] == 'C1_C2'].iloc[0]
            up_ee, down_es = map(int, c1c2_row['junction_coord'].split('-'))

            c1a_row = group[group['junction_name'] == 'C1_A'].iloc[0]
            exon_start, exon_end = map(int, c1a_row['spliced_with_coord'].split('-'))

            uniform_ID = f"{gene_id};SE:{chrom}:{up_ee}-{exon_start}:{exon_end}-{down_es}:{strand}"

        except (IndexError, ValueError):
            uniform_ID = pd.NA
        group = group.copy()
        group['uniform_ID'] = uniform_ID
        records.append(group)
    return pd.concat(records, ignore_index=True)


def parse_MAJIQ_A3SS_A5SS(df, event_type):
    """
    Parse MAJIQ A3SS and A5SS events to build uniform_ID.
    """
    records = []
    for event_id, group in df.groupby('event_id'):
        group = group.reset_index(drop=True)
        raw_gene_id = str(group.loc[0, 'gene_id'])
        match = re.search(r'ENSG\d+', raw_gene_id)
        gene_id = match.group(0) if match else raw_gene_id
        strand = group.loc[0, 'strand']
        chrom = group.loc[0, 'seqid'].replace('seqid', '')

        try:
            Proximal_row = group[group['junction_name'] == 'Proximal'].iloc[0]
            Pro_jun_coord_1, Pro_jun_coord_2 = [
                int(x) if x.isdigit() else x
                for x in Proximal_row['junction_coord'].split('-') if x.strip()
            ]
            Distal_row = group[group['junction_name'] == 'Distal'].iloc[0]
            Dis_jun_coord_1, Dis_jun_coord_2 = [
                int(x) if x.isdigit() else x
                for x in Distal_row['junction_coord'].split('-') if x.strip()
            ]

            if event_type == "A3SS":
                uniform_id = f"{gene_id};A3:{chrom}:{Pro_jun_coord_1}-{Pro_jun_coord_2}:{Dis_jun_coord_1}-{Dis_jun_coord_2}:{strand}"
            elif event_type == "A5SS":
                uniform_id = f"{gene_id};A5:{chrom}:{Pro_jun_coord_1}-{Pro_jun_coord_2}:{Dis_jun_coord_1}-{Dis_jun_coord_2}:{strand}"
        except (IndexError, ValueError):
            uniform_id = pd.NA
        group = group.copy()
        group['uniform_ID'] = uniform_id
        records.append(group)

    return pd.concat(records, ignore_index=True)


def parse_MAJIQ_AF_AL(df, event_type):
    """
    Parse MAJIQ AF (Alternative First Exon) and AL (Alternative Last Exon) events.
    """
    records = []
    for event_id, group in df.groupby('event_id'):
        group = group.reset_index(drop=True)
        raw_gene_id = str(group.loc[0, 'gene_id'])
        match = re.search(r'ENSG\d+', raw_gene_id)
        gene_id = match.group(0) if match else raw_gene_id
        strand = group.loc[0, 'strand']
        chrom = group.loc[0, 'seqid'].replace('seqid', '')

        try:
            Proximal_row = group[group['junction_name'] == 'Proximal'].iloc[0]
            Pro_jun_coord_1, Pro_jun_coord_2 = [
                int(x) if x.isdigit() else x
                for x in Proximal_row['junction_coord'].split('-') if x.strip()
            ]
            Pro_exon_start, Pro_exon_end = [
                int(x) if x.isdigit() else x
                for x in Proximal_row['spliced_with_coord'].split('-') if x.strip()
            ]

            Distal_row = group[group['junction_name'] == 'Distal'].iloc[0]
            Dis_jun_coord_1, Dis_jun_coord_2 = [
                int(x) if x.isdigit() else x
                for x in Distal_row['junction_coord'].split('-') if x.strip()
            ]
            Dis_exon_start, Dis_exon_end = [
                int(x) if x.isdigit() else x
                for x in Distal_row['spliced_with_coord'].split('-') if x.strip()
            ]

            if event_type == "AF":
                if strand == "+":
                    uniform_id = f"{gene_id};AF:{chrom}:{Dis_exon_start}:{Dis_exon_end}-{Pro_jun_coord_2}:{Pro_exon_start}:{Pro_exon_end}-{Dis_jun_coord_2}:{strand}"
                else:
                    uniform_id = f"{gene_id};AF:{chrom}:{Pro_jun_coord_1}-{Pro_exon_start}:{Pro_exon_end}:{Dis_jun_coord_1}-{Dis_exon_start}:{Dis_exon_end}:{strand}"
            elif event_type == "AL":
                if strand == "+":
                    uniform_id = f"{gene_id};AL:{chrom}:{Pro_jun_coord_1}-{Pro_exon_start}:{Pro_exon_end}:{Dis_jun_coord_1}-{Dis_exon_start}:{Dis_exon_end}:{strand}"
                else:
                    uniform_id = f"{gene_id};AL:{chrom}:{Dis_exon_start}:{Dis_exon_end}-{Pro_jun_coord_2}:{Pro_exon_start}:{Pro_exon_end}-{Dis_jun_coord_2}:{strand}"
        except (IndexError, ValueError):
            uniform_id = pd.NA
        group = group.copy()
        group['uniform_ID'] = uniform_id
        records.append(group)

    return pd.concat(records, ignore_index=True)


def parse_MAJIQ_RI(df, event_type):
    """
    Parse MAJIQ RI (Retained Intron) events.
    """
    records = []
    for event_id, group in df.groupby('event_id'):
        group = group.reset_index(drop=True)
        raw_gene_id = str(group.loc[0, 'gene_id'])
        match = re.search(r'ENSG\d+', raw_gene_id)
        gene_id = match.group(0) if match else raw_gene_id
        strand = group.loc[0, 'strand']
        chrom = group.loc[0, 'seqid'].replace('seqid', '')

        try:
            if group['spliced_with'].iloc[0] == "C1":
                intron_row = group[group['junction_name'] == 'C2_C1_intron'].iloc[0]
                intron_start, intron_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['junction_coord'].split('-') if x.strip()
                ]
                c1_start, c1_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['spliced_with_coord'].split('-') if x.strip()
                ]
                c2_start, c2_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['reference_exon_coord'].split('-') if x.strip()
                ]
                if strand == "+":
                    uniform_id = f"{gene_id};RI:{chrom}:{c1_start}:{intron_start-1}-{intron_end+1}:{c2_end}:{strand}"
                else:
                    uniform_id = f"{gene_id};RI:{chrom}:{c2_start}:{intron_start-1}-{intron_end+1}:{c1_end}:{strand}"
            elif group['spliced_with'].iloc[0] == "C2":
                intron_row = group[group['junction_name'] == 'C1_C2_intron'].iloc[0]
                intron_start, intron_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['junction_coord'].split('-') if x.strip()
                ]
                c2_start, c2_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['spliced_with_coord'].split('-') if x.strip()
                ]
                c1_start, c1_end = [
                    int(x) if x.isdigit() else x
                    for x in intron_row['reference_exon_coord'].split('-') if x.strip()
                ]
                if strand == "+":
                    uniform_id = f"{gene_id};RI:{chrom}:{c1_start}:{intron_start-1}-{intron_end+1}:{c2_end}:{strand}"
                else:
                    uniform_id = f"{gene_id};RI:{chrom}:{c2_start}:{intron_start-1}-{intron_end+1}:{c1_end}:{strand}"
        except (IndexError, ValueError):
            uniform_id = pd.NA
        group = group.copy()
        group['uniform_ID'] = uniform_id
        records.append(group)

    return pd.concat(records, ignore_index=True)


def parse_MAJIQ_MX(df, event_type):
    """
    Parse MAJIQ MX (Mutually Exclusive Exon) events.
    """
    records = []
    for event_id, group in df.groupby('event_id'):
        group = group.reset_index(drop=True)
        raw_gene_id = str(group.loc[0, 'gene_id'])
        match = re.search(r'ENSG\d+', raw_gene_id)
        gene_id = match.group(0) if match else raw_gene_id
        strand = group.loc[0, 'strand']
        chrom = str(group.loc[0, 'seqid']).replace('seqid', '')

        try:
            c1a1_row = group[group['junction_name'] == 'C1_A1'].iloc[0]
            coord_1, coord_2 = map(int, c1a1_row['junction_coord'].split('-'))

            c2a2_row = group[group['junction_name'] == 'C2_A2'].iloc[0]
            coord_3, coord_4 = map(int, c2a2_row['junction_coord'].split('-'))

            exon1_row = group[group['spliced_with'] == 'A1'].iloc[0]
            exon1_start, exon1_end = map(int, exon1_row['spliced_with_coord'].split('-'))

            exon2_row = group[group['spliced_with'] == 'A2'].iloc[0]
            exon2_start, exon2_end = map(int, exon2_row['spliced_with_coord'].split('-'))

            if strand == "+":
                uniform_ID = f"{gene_id};MX:{chrom}:{coord_1}-{exon1_start}:{exon1_end}-{coord_4}:{coord_1}-{exon2_start}:{exon2_end}-{coord_4}:{strand}"
            else:
                uniform_ID = f"{gene_id};MX:{chrom}:{coord_3}-{exon2_start}:{exon2_end}-{coord_2}:{coord_3}-{exon1_start}:{exon1_end}-{coord_2}:{strand}"
        except (IndexError, ValueError):
            uniform_ID = pd.NA
        group = group.copy()
        group['uniform_ID'] = uniform_ID
        records.append(group)
    return pd.concat(records, ignore_index=True)
