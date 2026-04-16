import numpy as np
import pandas as pd


def assign_event_id(
    df,
    coord_cols,
    group_cols=('gene', 'chrom', 'strand'),
    tool_col='tool',
    tol=1
):
    """
    Assign event_id by clustering splicing events using coordinate similarity
    and union-find algorithm. Ensures no duplicate tool within a cluster.
    """

    n = len(df)

    # ---------- parse coordinates ----------
    k = len(coord_cols)
    coords = np.empty((n, k * 2), dtype=np.int64)
    valid = np.ones(n, dtype=bool)

    for ci, c in enumerate(coord_cols):
        col = df[c].astype(str).values
        for i, v in enumerate(col):

            if 'na' in v.lower():
                valid[i] = False
                continue

            if '-' in v:
                s, e = v.split('-')
                coords[i, 2 * ci] = int(s)
                coords[i, 2 * ci + 1] = int(e)
            else:
                x = int(v)
                coords[i, 2 * ci] = x
                coords[i, 2 * ci + 1] = x

    tools = df[tool_col].values

    # ---------- union-find ----------
    parent = np.arange(n)

    # track tools per component
    comp_tools = [{tools[i]} for i in range(n)]

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra == rb:
            return

        # prevent same tool in one component
        if len(comp_tools[ra] & comp_tools[rb]) > 0:
            return

        parent[rb] = ra
        comp_tools[ra].update(comp_tools[rb])

    # ---------- grouping ----------
    groups = df.groupby(list(group_cols)).indices

    for gidx in groups.values():

        gidx = np.array([i for i in gidx if valid[i]])
        if len(gidx) == 0:
            continue

        order = np.argsort(coords[gidx, 0])
        idx = gidx[order]

        subcoords = coords[idx]

        m = len(idx)
        left = 0

        for i in range(m):

            while subcoords[i, 0] - subcoords[left, 0] > tol:
                left += 1

            for j in range(left, i):
                if np.all(np.abs(subcoords[i] - subcoords[j]) <= tol):
                    union(idx[i], idx[j])

    # ---------- compress ----------
    roots = np.array([find(i) for i in range(n)])
    _, event_id = np.unique(roots, return_inverse=True)

    return event_id


def fix_majiq_na(
    df,
    coord_cols,
    event_col='event_id',
    tool_col='tool',
    tol=1
):
    """
    Fix MAJIQ events with missing coordinates by matching them
    to nearest non-MAJIQ events.
    """

    df = df.copy()

    n = len(df)
    k = len(coord_cols)

    coords = np.zeros((n, k * 2), dtype=np.int64)
    mask = np.ones((n, k * 2), dtype=bool)

    for ci, c in enumerate(coord_cols):
        col = df[c].astype(str).values

        for i, v in enumerate(col):

            if '-' in v:
                s, e = v.split('-')

                if 'na' in s.lower():
                    mask[i, 2 * ci] = False
                else:
                    coords[i, 2 * ci] = int(s)

                if 'na' in e.lower():
                    mask[i, 2 * ci + 1] = False
                else:
                    coords[i, 2 * ci + 1] = int(e)

            else:
                if 'na' in v.lower():
                    mask[i, 2 * ci] = False
                    mask[i, 2 * ci + 1] = False
                else:
                    x = int(v)
                    coords[i, 2 * ci] = x
                    coords[i, 2 * ci + 1] = x

    tools = df[tool_col].values
    event = df[event_col].values

    tools_lower = np.char.lower(tools.astype(str))
    majiq_mask = np.char.startswith(tools_lower, "majiq")

    majiq_idx = np.where(majiq_mask & (~mask.all(axis=1)))[0]
    other_idx = np.where(~majiq_mask)[0]

    for i in majiq_idx:

        valid = mask[i]
        diff = np.abs(coords[other_idx][:, valid] - coords[i, valid])

        match = np.all(diff <= tol, axis=1)

        if np.any(match):
            j = other_idx[np.where(match)[0][0]]
            event[i] = event[j]

    df[event_col] = event
    return df

def merge_splicing_events(dfs, event_type):

    has_whippet = any(k.lower() == "whippet" for k in dfs)
    df = pd.concat(
        [d.assign(tool=k) for k, d in dfs.items()],
        ignore_index=True
    )
    all_whippet = all("whippet" in k.lower() for k in dfs)

    if all_whippet:
        df["event_id"] = pd.factorize(df["uniform_ID"])[0] + 1

    # 计算 frequency
        freq = df.groupby("event_id")["tool"].nunique()
        df["frequency"] = df["event_id"].map(freq)
        tmp = df["uniform_ID"].str.split(";", expand=True)
        p2 = tmp[1].str.split(",", expand=True)
        df["type"] = p2[0]
        df["gene"] = tmp[0]
        df["chrom"] = p2[1]
        df["strand"] = p2[2]
        return df
      
    tmp = df["uniform_ID"].str.split(";", expand=True)

    df["gene"] = tmp[0]

    p2 = tmp[1].str.split(",", expand=True)
    df["type"] = p2[0]
    df["chrom"] = p2[1]
    df["strand"] = p2[2]

    coords = tmp[2].str.split(":", expand=True)
    coord_cols = [f"c{i}" for i in range(coords.shape[1])]
    coords.columns = coord_cols

    df = pd.concat([df, coords], axis=1)

    if has_whippet:
        main = df[~df["tool"].str.lower().str.contains("whippet", na=False)].copy()
        whip = df[df["tool"].str.lower().str.contains("whippet", na=False)].copy()
    else:
        main = df.copy()

    # main event_id
    main['event_id']  = assign_event_id(main,coord_cols)
    main = fix_majiq_na(main, coord_cols)
    #main["event_id"] = pd.factorize(main["uniform_ID"])[0] + 1
    #main = propagate_event_id(main, coord_cols=coord_cols)

    if not has_whippet:

        freq = main.groupby("event_id")["tool"].nunique()
        main["frequency"] = main["event_id"].map(freq)

        return main

    whip_keys = []
    main_keys = []

    def build_key(d, coord):
        return (
            d["gene"] + "|" +
            d["type"] + "|" +
            d["chrom"] + "|" +
            d["strand"] + "|" +
            d[coord]
        )

    # SE / AF / AL
    if event_type in ["SE"]:

        whip_keys.append(pd.DataFrame({
            "uniform_ID":whip.uniform_ID,
            "key":build_key(whip,"c1")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(main,"c1"),
            "event_id":main.event_id
        }))

    # A5SS
    elif event_type=="A5SS":

        w = whip[(whip.strand=="+")]
        m = main[(main.strand=="+")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c0")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c0"),
            "event_id":m.event_id
        }))

        w = whip[(whip.strand=="-")]
        m = main[(main.strand=="-")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c1")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c1"),
            "event_id":m.event_id
        }))
  
    # A3SS
    elif event_type=="A3SS":

        w = whip[(whip.strand=="+")]
        m = main[(main.strand=="+")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c1")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c1"),
            "event_id":m.event_id
        }))

        w = whip[(whip.strand=="-")]
        m = main[(main.strand=="-")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c0")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c0"),
            "event_id":m.event_id
        }))
        
    elif event_type=="AF":

        w = whip[(whip.strand=="+")]
        m = main[(main.strand=="+")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c2")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c2"),
            "event_id":m.event_id
        }))

        w = whip[(whip.strand=="-")]
        m = main[(main.strand=="-")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c1")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c1"),
            "event_id":m.event_id
        }))
        
    elif event_type=="AL":

        w = whip[(whip.strand=="+")]
        m = main[(main.strand=="+")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c1")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c1"),
            "event_id":m.event_id
        }))

        w = whip[(whip.strand=="-")]
        m = main[(main.strand=="-")]

        whip_keys.append(pd.DataFrame({
            "uniform_ID":w.uniform_ID,
            "key":build_key(w,"c2")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(m,"c2"),
            "event_id":m.event_id
        }))
    # RI
    elif event_type=="RI":
        whip["c"] = whip["c1"] + "-" + whip["c3"]
        main["c"] = main["c1"] + "-" + main["c3"]
        whip_keys.append(pd.DataFrame({
            "uniform_ID":whip.uniform_ID,
            "key":build_key(whip,"c")
        }))

        main_keys.append(pd.DataFrame({
            "key":build_key(main,"c"),
            "event_id":main.event_id
        }))

    whip_keys = pd.concat(whip_keys)
    main_keys = pd.concat(main_keys)

    match = whip_keys.merge(main_keys,on="key",how="left")
    match = match.drop_duplicates(["uniform_ID","event_id"])
    whip_match = whip.merge(
        match[["uniform_ID","event_id"]],
        on="uniform_ID",
        how="left"
    )

    max_id = main.event_id.max()
    mask = whip_match.event_id.isna()

    whip_match.loc[mask,"event_id"] = range(
        max_id+1,
        max_id+1+mask.sum()
    )

    df_final = pd.concat([main,whip_match],ignore_index=True)

    freq = df_final.groupby("event_id")["tool"].nunique()
    df_final["frequency"] = df_final["event_id"].map(freq)

    return df_final
