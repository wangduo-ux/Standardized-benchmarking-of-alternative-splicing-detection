import pandas as pd


# ===== SE =====
def parse_rMATS_SE(df, event_type="SE"):
    df = df.copy()

    gene_id = df["GeneID"].astype(str).str.strip('"')
    chrom = df["chr"].astype(str).str.replace(r"^chr", "", regex=True)
    strand = df["strand"].astype(str)

    # safer numeric conversion (column-wise, not all-or-nothing)
    up_ee = pd.to_numeric(df["upstreamEE"], errors='coerce').astype('Int64')
    exon_start = pd.to_numeric(df["exonStart_0base"], errors='coerce').astype('Int64') + 1
    exon_end = pd.to_numeric(df["exonEnd"], errors='coerce').astype('Int64')
    down_es = pd.to_numeric(df["downstreamES"], errors='coerce').astype('Int64') + 1

    df["uniform_ID"] = (
        gene_id + ";SE," + chrom + "," + strand + ";" +
        up_ee.astype(str) + ":" +
        exon_start.astype(str) + "-" + exon_end.astype(str) + ":" +
        down_es.astype(str)
    )

    return df


# ===== A3SS / A5SS =====
def parse_rMATS_A3SS_A5SS(df, event_type):
    if event_type not in ["A3SS", "A5SS"]:
        raise ValueError(f"Invalid event_type for rMATS: {event_type}")

    def generate_ID(row):
        try:
            gene_id = str(row["GeneID"]).strip('"')
            chrom = str(row["chr"])
            strand = row["strand"]

            if chrom.lower().startswith("chr"):
                chrom = chrom[3:]

            longExonStart = int(row["longExonStart_0base"])
            longExonEnd = int(row["longExonEnd"])
            shortES = int(row["shortES"])
            shortEE = int(row["shortEE"])
            flankingES = int(row["flankingES"])
            flankingEE = int(row["flankingEE"])

            if event_type == "A3SS":
                if strand == "+":
                    return f"{gene_id};A3SS,{chrom},{strand};{flankingEE}:{longExonStart+1}-{shortES+1}"
                else:
                    return f"{gene_id};A3SS,{chrom},{strand};{shortEE}-{longExonEnd}:{flankingES+1}"

            elif event_type == "A5SS":
                if strand == "+":
                    return f"{gene_id};A5SS,{chrom},{strand};{shortEE}-{longExonEnd}:{flankingES+1}"
                else:
                    return f"{gene_id};A5SS,{chrom},{strand};{flankingEE}:{longExonStart+1}-{shortES+1}"

        except Exception as e:
            raise ValueError(f"rMATS {event_type} parse error: {row.to_dict()} | {e}")

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df


# ===== RI =====
def parse_rMATS_RI(df, event_type="RI"):

    def generate_ID(row):
        try:
            gene_id = str(row["GeneID"]).strip('"')
            chrom = str(row["chr"])
            strand = row["strand"]

            if chrom.lower().startswith("chr"):
                chrom = chrom[3:]

            upstreamES = int(row["upstreamES"])
            upstreamEE = int(row["upstreamEE"])
            downstreamES = int(row["downstreamES"])
            downstreamEE = int(row["downstreamEE"])

            return f"{gene_id};RI,{chrom},{strand};{upstreamES+1}:{upstreamEE}-{downstreamES+1}:{downstreamEE}:{upstreamEE}-{downstreamES+1}"

        except Exception as e:
            raise ValueError(f"rMATS RI parse error: {row.to_dict()} | {e}")

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df


# ===== MX =====
def parse_rMATS_MX(df, event_type="MX"):

    def generate_ID(row):
        try:
            gene_id = str(row["GeneID"]).strip('"')
            chrom = str(row["chr"])
            strand = row["strand"]

            if chrom.lower().startswith("chr"):
                chrom = chrom[3:]

            exon1_start = int(row["1stExonStart_0base"])
            exon1_end = int(row["1stExonEnd"])
            exon2_start = int(row["2ndExonStart_0base"])
            exon2_end = int(row["2ndExonEnd"])
            upstreamEE = int(row["upstreamEE"])
            downstreamES = int(row["downstreamES"])

            return (
                f"{gene_id};MX,{chrom},{strand};"
                f"{upstreamEE}:{exon1_start+1}-{exon1_end}:"
                f"{exon2_start+1}-{exon2_end}:{downstreamES+1}"
            )

        except Exception as e:
            raise ValueError(f"rMATS MX parse error: {row.to_dict()} | {e}")

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df
