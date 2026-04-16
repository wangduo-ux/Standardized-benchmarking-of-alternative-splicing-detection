import pandas as pd


def extract_dpsi(df, sw, ev):
    """
    Extract delta PSI / PSI change values from different splicing tools
    and normalize into a unified 'value' column.
    """

    if 'suppa2' in sw.lower():

        df = df.loc[:, ['uniform_ID'] + [col for col in df.columns if 'dPSI' in col]]
        df['value'] = df.filter(like='dPSI')

    elif 'rmats' in sw.lower():

        df = df[['uniform_ID', 'IncLevelDifference']]
        df = df.rename(columns={'IncLevelDifference': 'value'})

    elif 'psi-sigma' in sw.lower():

        if ev == "MX":
            df = df.loc[df.groupby('uniform_ID')['Target Exon'].idxmin()]

        df = df[['uniform_ID', 'ΔPSI (%)']]
        df = df.rename(columns={'ΔPSI (%)': 'value'})
        df['value'] = df['value'] / 100

    elif 'majiq' in sw.lower():

        if ev == "SE":
            df = df[df['spliced_with'] == 'A']

        elif ev in ["A3SS", "A5SS"]:
            df = df[df['junction_name'] == 'Proximal']

        elif ev in ["AF", "AL"]:
            df = df[df['junction_name'] == 'Distal']

        elif ev == "MX":

            pos = df[(df["strand"] == "+") & (df["junction_name"].isin(["C1_A1"]))]
            neg = df[(df["strand"] == "-") & (df["junction_name"].isin(["C1_A2"]))]

            df = pd.concat([pos, neg], ignore_index=True)

        elif ev == "RI":
            df = df[df['junction_name'].str.contains('intron')]

        dpsi_col = next(col for col in df.columns if 'dpsi' in col)
        df = df[['uniform_ID', dpsi_col]]

        df = (
            df.groupby('uniform_ID')[dpsi_col]
            .mean()
            .reset_index()
        )

        df = df.rename(columns={dpsi_col: 'value'})

    elif 'spladder' in sw.lower():

        df = df[['uniform_ID', 'dPSI']]
        df = df.rename(columns={'dPSI': 'value'})

    elif 'whippet' in sw.lower():

        df = df[['uniform_ID', 'DeltaPsi']]
        df = df.rename(columns={'DeltaPsi': 'value'})

        if ev in ["AF", "AL"]:
            df["value"] = -df["value"]

    else:
        raise ValueError(f"Unsupported software: {sw}")

    # ---------- final cleanup ----------
    df = df[['uniform_ID', 'value']]
    df = df.dropna(subset=["value"])

    return df
