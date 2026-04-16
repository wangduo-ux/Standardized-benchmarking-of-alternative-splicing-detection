import pandas as pd


# ---------- helper ----------

def clean_chrom(chrom):
    chrom = str(chrom)
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


# ---------- SE ----------

def parse_Spladder_SE(df, event_type="SE"):
    """
    Parse Spladder SE events and generate uniform_ID.
    """

    def generate_ID(row):
        gene_id = row["gene_id"]
        chrom = clean_chrom(row["chrm"])
        strand = row["strand"]

        try:
            up_es, up_ee, start, end, down_es, down_ee = \
                map(int, row["exon_pos"].replace(":", "-").split("-"))
        except Exception:
            return pd.NA

        start += 1
        down_es += 1

        return f"{gene_id};SE,{chrom},{strand};{up_ee}:{start}-{end}:{down_es}"

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df


# ---------- A3SS / A5SS ----------

def parse_Spladder_A3SS_A5SS(df, event_type):
    """
    Parse Spladder A3SS / A5SS events.
    """

    def generate_ID(row):
        gene_id = row["gene_id"]
        chrom = clean_chrom(row["chrm"])
        strand = row["strand"]

        try:
            exon1_es, exon1_ee, exon2_es, exon2_ee, exon3_es, exon3_ee = \
                map(int, row["exon_pos"].replace(":", "-").split("-"))
        except Exception:
            return pd.NA

        if event_type == "A3SS":

            if strand == "+":
                return f"{gene_id};A3SS,{chrom},{strand};{exon1_ee}:{exon2_es+1}-{exon3_es+1}"
            else:
                return f"{gene_id};A3SS,{chrom},{strand};{exon1_ee}-{exon2_ee}:{exon3_es+1}"

        elif event_type == "A5SS":

            if strand == "+":
                return f"{gene_id};A5SS,{chrom},{strand};{exon1_ee}-{exon2_ee}:{exon3_es+1}"
            else:
                return f"{gene_id};A5SS,{chrom},{strand};{exon1_ee}:{exon2_es+1}-{exon3_es+1}"

        return pd.NA

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df


# ---------- RI ----------

def parse_Spladder_RI(df, event_type="RI"):
    """
    Parse Spladder RI events.
    """

    def generate_ID(row):
        gene_id = row["gene_id"]
        chrom = clean_chrom(row["chrm"])
        strand = row["strand"]

        try:
            exon1_es, exon1_ee, exon2_es, exon2_ee, exon3_es, exon3_ee = \
                map(int, row["exon_pos"].replace(":", "-").split("-"))
        except Exception:
            return pd.NA

        return (
            f"{gene_id};RI,{chrom},{strand};"
            f"{exon1_es+1}:{exon1_ee}-{exon3_es+1}:{exon3_ee}:{exon1_ee}-{exon3_es+1}"
        )

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df


# ---------- MX ----------

def parse_Spladder_MX(df, event_type="MX"):
    """
    Parse Spladder MX events.
    """

    def generate_ID(row):
        gene_id = row["gene_id"]
        chrom = clean_chrom(row["chrm"])
        strand = row["strand"]

        try:
            up_es, up_ee, e1s, e1e, e2s, e2e, down_es, down_ee = \
                map(int, row["exon_pos"].replace(":", "-").split("-"))
        except Exception:
            return pd.NA

        e1s += 1
        e2s += 1
        down_es += 1

        return f"{gene_id};MX,{chrom},{strand};{up_ee}:{e1s}-{e1e}:{e2s}-{e2e}:{down_es}"

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df
