import pandas as pd
import re


# ---------- helper functions ----------

def clean_chrom(chrom):
    chrom = str(chrom)
    if chrom.lower().startswith("chr"):
        return chrom[3:]
    return chrom


def safe_gene_id(x):
    x = str(x)
    return x.split(":", 1)[-1].strip('"')


def safe_extract(series, pattern):
    out = series.str.extract(pattern)
    return out.fillna("na")


def parse_coord(x):
    m = re.match(r'(-?\d+|na)-(-?\d+|na)', str(x))
    if not m:
        return "na", "na"
    return m.group(1), m.group(2)


# ---------- SE ----------

def parse_MAJIQ_SE(df, event_type="SE"):

    df = df.copy()

    c1c2 = df[df["junction_name"] == "C1_C2"][["event_id", "junction_coord"]].copy()
    c1c2[["up_ee", "down_es"]] = safe_extract(
        c1c2["junction_coord"], r'(-?\d+|na)-(-?\d+|na)'
    )

    c1a = df[df["junction_name"] == "C1_A"][["gene_id", "strand", "event_id", "seqid", "spliced_with_coord"]].copy()
    c1a[["exon_start", "exon_end"]] = safe_extract(
        c1a["spliced_with_coord"], r'(-?\d+|na)-(-?\d+|na)'
    )

    c1a[["exon_start", "exon_end"]] = c1a[["exon_start", "exon_end"]].replace("-1", "na")

    c1a["gene_id"] = c1a["gene_id"].apply(safe_gene_id)
    c1a["strand"] = c1a["strand"].str[0]
    c1a["seqid"] = c1a["seqid"].apply(clean_chrom)

    coords = c1c2.merge(c1a, on="event_id", how="inner")

    coords["uniform_ID"] = (
        coords["gene_id"] + ";SE," + coords["seqid"] + "," + coords["strand"] + ";" +
        coords["up_ee"] + ":" +
        coords["exon_start"] + "-" + coords["exon_end"] + ":" +
        coords["down_es"]
    )

    df = df.merge(coords[["event_id", "uniform_ID"]], on="event_id", how="left")

    return df


# ---------- A3SS / A5SS ----------

def parse_MAJIQ_A3SS_A5SS(df, event_type):

    df = df.copy()

    proximal = df[df["junction_name"] == "Proximal"][["gene_id", "strand", "event_id", "seqid", "junction_coord"]].copy()
    distal = df[df["junction_name"] == "Distal"][["event_id", "junction_coord"]].groupby("event_id", as_index=False).first()

    proximal["gene_id"] = proximal["gene_id"].apply(safe_gene_id)
    proximal["strand"] = proximal["strand"].str[0]
    proximal["seqid"] = proximal["seqid"].apply(clean_chrom)

    proximal[["Pro1", "Pro2"]] = safe_extract(proximal["junction_coord"], r'(-?\d+|na)-(-?\d+|na)')
    distal[["Dis1", "Dis2"]] = safe_extract(distal["junction_coord"], r'(-?\d+|na)-(-?\d+|na)')

    coords = proximal.merge(distal, on="event_id", how="inner")

    def build_uniform_id(row):
        gene_id = row["gene_id"]
        strand = row["strand"]
        chrom = row["seqid"]
        Pro1, Pro2, Dis1, Dis2 = row["Pro1"], row["Pro2"], row["Dis1"], row["Dis2"]

        if event_type == "A3SS":
            if strand == "+":
                return f"{gene_id};A3SS,{chrom},{strand};{Pro1}:{Pro2}-{Dis2}"
            else:
                return f"{gene_id};A3SS,{chrom},{strand};{Dis1}-{Pro1}:{Pro2}"
        elif event_type == "A5SS":
            if strand == "+":
                return f"{gene_id};A5SS,{chrom},{strand};{Dis1}-{Pro1}:{Pro2}"
            else:
                return f"{gene_id};A5SS,{chrom},{strand};{Pro1}:{Pro2}-{Dis2}"

    coords["uniform_ID"] = coords.apply(build_uniform_id, axis=1)

    df = df.merge(coords[["event_id", "uniform_ID"]], on="event_id", how="left")

    return df


# ---------- AF / AL ----------

def parse_MAJIQ_AF_AL(df, event_type):

    df = df.copy()

    proximal = df[df["junction_name"] == "Proximal"][["gene_id", "strand", "event_id", "seqid", "junction_coord", "spliced_with_coord"]].copy()
    distal = df[df["junction_name"] == "Distal"][["event_id", "junction_coord", "spliced_with_coord"]].copy()

    proximal["gene_id"] = proximal["gene_id"].apply(safe_gene_id)
    proximal["seqid"] = proximal["seqid"].apply(clean_chrom)

    proximal[["Pro_j1", "Pro_j2"]] = safe_extract(proximal["junction_coord"], r'(-?\d+|na)-(-?\d+|na)')
    proximal[["Pro_ex_start", "Pro_ex_end"]] = safe_extract(proximal["spliced_with_coord"], r'(-?\d+|na)-(-?\d+|na)')

    distal[["Dis_j1", "Dis_j2"]] = safe_extract(distal["junction_coord"], r'(-?\d+|na)-(-?\d+|na)')
    distal[["Dis_ex_start", "Dis_ex_end"]] = safe_extract(distal["spliced_with_coord"], r'(-?\d+|na)-(-?\d+|na)')

    coords = proximal.merge(distal, on="event_id", how="inner")

    def build_uniform_id(row):
        gene_id = row["gene_id"]
        strand = row["strand"]
        chrom = row["seqid"]

        Pro_j1, Pro_j2 = row["Pro_j1"], row["Pro_j2"]
        Pro_ex_start, Pro_ex_end = row["Pro_ex_start"], row["Pro_ex_end"]
        Dis_j1, Dis_j2 = row["Dis_j1"], row["Dis_j2"]
        Dis_ex_start, Dis_ex_end = row["Dis_ex_start"], row["Dis_ex_end"]

        if event_type == "AF":
            if strand == "+":
                return f"{gene_id};AF,{chrom},{strand};{Dis_ex_start}-{Dis_j1}:{Dis_j2}:{Pro_ex_start}-{Pro_j1}:{Pro_j2}"
            else:
                return f"{gene_id};AF,{chrom},{strand};{Pro_j1}:{Pro_j2}-{Pro_ex_end}:{Dis_j1}:{Dis_j2}-{Dis_ex_end}"

        elif event_type == "AL":
            if strand == "+":
                return f"{gene_id};AL,{chrom},{strand};{Pro_j1}:{Pro_j2}-{Pro_ex_end}:{Dis_j1}:{Dis_j2}-{Dis_ex_end}"
            else:
                return f"{gene_id};AL,{chrom},{strand};{Dis_ex_start}-{Dis_j1}:{Dis_j2}:{Pro_ex_start}-{Pro_j1}:{Pro_j2}"

    coords["uniform_ID"] = coords.apply(build_uniform_id, axis=1)

    df = df.merge(coords[["event_id", "uniform_ID"]], on="event_id", how="left")

    return df


# ---------- RI ----------

def parse_MAJIQ_RI(df, event_type="RI"):
    
    intron = df[df["junction_name"].isin(["C2_C1_intron", "C1_C2_intron"])][["gene_id","strand","event_id","seqid","junction_coord","reference_exon_coord","spliced_with_coord","junction_name"]].copy()
    intron["gene_id"] = intron["gene_id"].str.replace(r".*:", "", regex=True)
    spliced = df[df["junction_name"].isin(["C2_C1_spliced", "C1_C2_spliced"])][["event_id","junction_coord","reference_exon_coord","spliced_with_coord","junction_name"]].copy()
    
    
    coords = intron.merge(spliced, on="event_id", how="inner",suffixes=("_intron", "_spliced"))
    def parse_coord(x):
        m = re.match(r'(-?\d+|na)-(-?\d+|na)', str(x))
        if not m:
            return None, None
        return m.group(1), m.group(2)
    # 构建 uniform_ID
    def build_uniform_id(row,event_type):
        gene_id = row["gene_id"]
        strand = row["strand"]
        chrom = str(row["seqid"]).replace("seqid", "")
        intron_start, intron_end = parse_coord(row["junction_coord_intron"])

        if row["junction_name_intron"] == "C2_C1_intron":
            c2_start, c2_end = parse_coord(row["reference_exon_coord_intron"])
            c1_start, c1_end = parse_coord(row["spliced_with_coord_intron"])
        else:
            c1_start, c1_end = parse_coord(row["reference_exon_coord_intron"])
            c2_start, c2_end = parse_coord(row["spliced_with_coord_intron"])

        c1_start, c1_end, c2_start, c2_end = [
            'na' if x == "-1" else x
            for x in [c1_start, c1_end, c2_start, c2_end]
        ]

        spliced_start, spliced_end = parse_coord(row["junction_coord_spliced"])

        if strand == "+":
            return f"{gene_id};RI,{chrom},{strand};{c1_start}:{intron_start}-{intron_end}:{c2_end}:{spliced_start}-{spliced_end}"
        else:
            return  f"{gene_id};RI,{chrom},{strand};{c2_start}:{intron_start}-{intron_end}:{c1_end}:{spliced_start}-{spliced_end}"
    
    
    coords["uniform_ID"] = coords.apply(lambda row: build_uniform_id(row, event_type), axis=1)
    coords = coords[["uniform_ID","event_id"]]
    # 合并回原 df
    df = df.merge(coords[["event_id","uniform_ID"]], on="event_id", how="left")
    return df
# ---------- MX ----------

def parse_MAJIQ_MX(df, event_type="MX"):

    df = df.copy()
    uniform_ids = pd.Series(index=df.index, dtype=object)

    for gid, g in df.groupby("event_id", sort=False):

        g = g.copy()
        g[["start", "end"]] = g["spliced_with_coord"].str.extract(r'(-?\d+|na)-(-?\d+|na)')
        g["start"] = g["start"].replace("na", None).astype(float)

        g_sorted = g.sort_values("start")

        try:
            gene_id = safe_gene_id(g["gene_id"].iloc[0])
            seqid = clean_chrom(g["seqid"].iloc[0])
            strand = g["strand"].iloc[0]

            first_type = g_sorted.iloc[0]["spliced_with"]

            if first_type == "A1":
                a1 = g[g.spliced_with == "A1"].iloc[0]
                a2 = g[g.spliced_with == "A2"].iloc[0]
                j1 = g[g.junction_name.str.contains("C1")].iloc[0]
                j2 = g[g.junction_name.str.contains("C2")].iloc[0]
            else:
                a1 = g[g.spliced_with == "A2"].iloc[0]
                a2 = g[g.spliced_with == "A1"].iloc[0]
                j1 = g[g.junction_name.str.contains("C1")].iloc[0]
                j2 = g[g.junction_name.str.contains("C2")].iloc[0]

        except Exception:
            continue

        a1_start, a1_end = parse_coord(a1.spliced_with_coord)
        a2_start, a2_end = parse_coord(a2.spliced_with_coord)
        jun1, jun2 = parse_coord(j1.junction_coord)
        jun3, jun4 = parse_coord(j2.junction_coord)

        if strand == "+":
            uid = f"{gene_id};MX,{seqid},{strand};{jun1}:{a1_start}-{a1_end}:{a2_start}-{a2_end}:{jun4}"
        else:
            uid = f"{gene_id};MX,{seqid},{strand};{jun3}:{a1_start}-{a1_end}:{a2_start}-{a2_end}:{jun2}"

        uniform_ids.loc[g.index] = uid

    df["uniform_ID"] = uniform_ids

    return df
