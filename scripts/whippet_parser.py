import pandas as pd


# ---------- helper ----------

def clean_chrom(chrom):
    chrom = str(chrom)
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


# ---------- main parser ----------

def parse_Whippet(df, event_type):
    """
    Parse Whippet output and generate uniform_ID.
    """

    def generate_ID(row):

        try:
            gene = str(row["Gene"])
            coord = str(row["Coord"])
            strand = row["Strand"]

            chrom, pos = coord.split(":")
            chrom = clean_chrom(chrom)

            start, end = map(int, pos.split("-"))

        except Exception:
            return pd.NA

        if event_type == "SE":
            return f"{gene};SE,{chrom},{strand};na:{start}-{end}:na"

        elif event_type == "AF":
            if strand == "+":
                return f"{gene};AF,{chrom},{strand};na-na:na:{start}-{end}:na"
            else:
                return f"{gene};AF,{chrom},{strand};na:{start}-{end}:na:na-na"

        elif event_type == "AL":
            if strand == "-":
                return f"{gene};AL,{chrom},{strand};na-na:na:{start}-{end}:na"
            else:
                return f"{gene};AL,{chrom},{strand};na:{start}-{end}:na:na-na"

        elif event_type == "A5SS":
            if strand == "+":
                return f"{gene};A5SS,{chrom},{strand};{start-1}-{end}:na"
            else:
                return f"{gene};A5SS,{chrom},{strand};na:{start}-{end+1}"

        elif event_type == "A3SS":
            if strand == "-":
                return f"{gene};A3SS,{chrom},{strand};{start-1}-{end}:na"
            else:
                return f"{gene};A3SS,{chrom},{strand};na:{start}-{end+1}"

        elif event_type == "RI":
            return f"{gene};RI,{chrom},{strand};na:{start-1}-{end+1}:na:{start-1}-{end+1}"

        return pd.NA

    df = df.copy()
    df["uniform_ID"] = df.apply(generate_ID, axis=1)

    return df
