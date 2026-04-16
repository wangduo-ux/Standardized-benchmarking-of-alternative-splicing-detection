import re
import pandas as pd
import itertools


    
    
def parse_psisigma_SE_A3SS_A5SS(df, event_type):

    def generate_ID(row):
        try:
            gene_id_raw = str(row["gene_id"])
            m = re.match(r'^"(ENSG[^"]+)"$', gene_id_raw)
            gene_id = m.group(1) if m else gene_id_raw.strip('"')

            Target_Exon = row["Target Exon"]
            Event_Region = row["Event Region"]
            Event_Type = row["Event Type"]
            strand = row["strand"]

            # --- Target exon ---
            m1 = re.match(r'^([^:]+):(\d+)-(\d+)$', Target_Exon)
            if not m1:
                return pd.NA

            chrom = m1.group(1)
            if chrom.lower().startswith('chr'):
                chrom = chrom[3:]

            exon_start = int(m1.group(2))
            exon_end = int(m1.group(3))

            # --- Event region ---
            m2 = re.match(r'^([^:]+):(\d+)-(\d+)$', Event_Region)
            if not m2:
                return pd.NA

            up_ee = int(m2.group(2))
            down_es = int(m2.group(3))

            # --- ID construction ---
            if Event_Type == "SES":
                return f"{gene_id};SE,{chrom},{strand};{up_ee-1}:{exon_start}-{exon_end}:{down_es+1}"

            elif Event_Type in ["A3SS", "TSS|A3SS"]:
                if strand == "+":
                    return f"{gene_id};A3SS,{chrom},{strand};{up_ee-1}:{exon_start}-{down_es+1}"
                else:
                    return f"{gene_id};A3SS,{chrom},{strand};{up_ee-1}-{exon_end}:{down_es+1}"

            elif Event_Type in ["A5SS", "TSS|A5SS"]:
                if strand == "+":
                    return f"{gene_id};A5SS,{chrom},{strand};{up_ee-1}-{exon_end}:{down_es+1}"
                else:
                    return f"{gene_id};A5SS,{chrom},{strand};{up_ee-1}:{exon_start}-{down_es+1}"

            return pd.NA

        except Exception as e:
            raise ValueError(f"PSI-Sigma parse error: {row.to_dict()} | {e}")

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df
  
  
def parse_psisigma_MX(df, event_type):

    df = df[df['Event Region'].map(df['Event Region'].value_counts()) > 1].copy()
    if df.empty:
        return df

    merged_dfs = []

    for event_region, group in df.groupby("Event Region"):
        try:
            chrom, coord = event_region.split(":")
            region_start, region_end = map(int, coord.split('-'))
        except Exception:
            continue

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

        for (i1, s1, e1), (i2, s2, e2) in itertools.combinations(exon_info, 2):
            if e1 < s2 or e2 < s1:

                if s1 < s2:
                    left = f"{s1}-{e1}"
                    right = f"{s2}-{e2}"
                    keep_index = i1
                else:
                    left = f"{s2}-{e2}"   # ✔ 修复 ":" bug
                    right = f"{s1}-{e1}"
                    keep_index = i2

                uniform_id = (
                    f"{gene_id};MX,{chrom},{strand};"
                    f"{region_start-1}:{left}:{right}:{region_end+1}"
                )

                merged_row = pd.concat([df.loc[i1], df.loc[i2]], axis=1).T
                merged_row["uniform_ID"] = uniform_id
                merged_row = merged_row.reset_index(drop=True)

                merged_dfs.append(merged_row)

    if not merged_dfs:
        return pd.DataFrame()

    return pd.concat(merged_dfs, ignore_index=True)
  
def parse_psisigma_RI(df, event_type):

    def generate_ID(row):
        try:
            gene_id_raw = str(row["gene_id"])
            m = re.match(r'^"(ENSG[^"]+)"$', gene_id_raw)
            gene_id = m.group(1) if m else gene_id_raw.strip('"')

            Target_Exon = row["Target Exon"]
            Event_Region = row["Event Region"]
            strand = row["strand"]

            m1 = re.match(r'^([^:]+):(\d+)-(\d+)$', Target_Exon)
            m2 = re.match(r'^([^:]+):(\d+)-(\d+)$', Event_Region)

            if not m1 or not m2:
                return pd.NA

            chrom = m1.group(1)
            if chrom.lower().startswith('chr'):
                chrom = chrom[3:]

            exon_start = int(m1.group(2))
            exon_end = int(m1.group(3))
            up_ee = int(m2.group(2))
            down_es = int(m2.group(3))

            return (
                f"{gene_id};RI,{chrom},{strand};"
                f"{up_ee}:{exon_start-1}-{exon_end+1}:{down_es}:{exon_start-1}-{exon_end+1}"
            )

        except Exception as e:
            raise ValueError(f"PSI-Sigma RI parse error: {row.to_dict()} | {e}")

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df
