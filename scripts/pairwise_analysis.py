import pandas as pd
import numpy as np
from itertools import combinations
import os


# =========================================================
# 1. Pairwise construction
# =========================================================

def build_pairwise_df(df):
    """
    Build pairwise PSI/dPSI comparison table across tools.
    Input must contain:
        ['event_id','uniform_ID','value','tool','type']
    """

    out = []

    for tp, sub in df.groupby("type"):

        tools = sub["tool"].unique()

        for t1, t2 in combinations(tools, 2):

            df1 = sub[sub.tool == t1][["event_id", "uniform_ID", "value"]].copy()
            df2 = sub[sub.tool == t2][["event_id", "uniform_ID", "value"]].copy()

            merged = df1.merge(
                df2,
                on="event_id",
                suffixes=(f"_{t1}", f"_{t2}")
            ).dropna()

            tmp = pd.DataFrame({
                "value1": merged[f"value_{t1}"],
                "value2": merged[f"value_{t2}"],
                "tool1": t1,
                "tool2": t2,
                "type": tp,
                "ID1": merged[f"uniform_ID_{t1}"],
                "ID2": merged[f"uniform_ID_{t2}"]
            })

            out.append(tmp)

    if len(out) == 0:
        return pd.DataFrame()

    return pd.concat(out, ignore_index=True)


# =========================================================
# 2. Pairwise correlation
# =========================================================

def compute_pairwise_pearson(df, output_dir=None, sample_name=None):
    """
    Compute Pearson/Spearman correlation for each tool pair.
    """

    results = []

    for tp, sub in df.groupby("type"):

        sub = sub.copy()
        sub["value1"] = pd.to_numeric(sub["value1"], errors="coerce")
        sub["value2"] = pd.to_numeric(sub["value2"], errors="coerce")
        sub = sub.dropna(subset=["value1", "value2"])

        for (t1, t2), g in sub.groupby(["tool1", "tool2"]):

            n = len(g)

            if n < 2:
                pearson_r = np.nan
                spearman_r = np.nan
            else:
                pearson_r = g["value1"].corr(g["value2"], method="pearson")
                spearman_r = g["value1"].corr(g["value2"], method="spearman")

            results.append({
                "type": tp,
                "tool1": t1,
                "tool2": t2,
                "pearson_r": pearson_r,
                "spearman_r": spearman_r,
                "n": n
            })

    results = pd.DataFrame(results)

    # optional IO
    if output_dir is not None and sample_name is not None:
        os.makedirs(output_dir, exist_ok=True)
        out_file = os.path.join(output_dir, f"{sample_name}.dPSI.correlation.txt")
        results.to_csv(out_file, sep="\t", index=False)

    return results


# =========================================================
# 3. Whippet overlap summary
# =========================================================


def summarize_whippet_overlap(df):
    if "Whippet_ID" not in df.columns:
        return None

    results = []

    id_cols = [c for c in df.columns if c.endswith("_ID") and c != "Whippet_ID"]

    for tp, g in df[df["Whippet_ID"].notna()].groupby("event_type"):

        total = g["Whippet_ID"].nunique()

        for c in id_cols:
            results.append({
                "event_type": tp,
                "software": c.replace("_ID",""),
                "overlap_whippet_n": g.loc[g[c].notna(),"Whippet_ID"].nunique(),
                "overlap_event_n": g.loc[g[c].notna(),"event_id"].nunique(),
                "total_whippet_n": total
            })

    return pd.DataFrame(results)
  


def build_upset_df(df):
    return (
        pd.concat([
            df.assign(present=1)
            .pivot_table(index="event_id", columns="tool", values="present", fill_value=0),
            df
            .pivot_table(index="event_id", columns="tool", values="uniform_ID", aggfunc="first")
            .add_suffix("_ID")
        ], axis=1)
        .reset_index()
    )
