"""
dpsi_plotting.py

Module functions:
1. extract_dpsi: Extract ΔPSI or similar values from alternative splicing data across different software.
2. plot_corr_dpsi: Plot correlation of ΔPSI values from multiple software using a PairGrid.
3. plot_all_events_grid: Combine correlation plots of multiple events into a grid and save as an image.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


def extract_dpsi(df, sw, ev):
    """
    Extract ΔPSI (or equivalent metric) and unify column name to 'value'.
    
    Parameters:
    - df: pd.DataFrame, input event data.
    - sw: str, software name ('SUPPA2', 'rMATS', 'MAJIQ', 'PSI-Sigma').
    - ev: str, event type ('SE', 'A3SS', 'A5SS', 'AF', 'AL', 'RI', 'MX').

    Returns:
    - pd.DataFrame: containing 'uniform_ID' and 'value' columns.
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
            pos = df[(df["strand"] == "+") & (df["junction_name"].isin(["C1_A1", "C2_A2"]))]
            neg = df[(df["strand"] == "-") & (df["junction_name"].isin(["C1_A2", "C2_A1"]))]
            df = pd.concat([pos, neg], ignore_index=True)
        elif ev == "RI":
            df = df[df['junction_name'].str.contains('intron')]
        dpsi_cols = [col for col in df.columns if 'dpsi' in col]
        df = df[['uniform_ID'] + dpsi_cols]
        df = df.groupby('uniform_ID')[dpsi_cols].mean().reset_index()
        df['value'] = df[dpsi_cols].mean(axis=1)
    return df


def plot_corr_dpsi(software_dfs):
    """
    Plot correlation of ΔPSI values from multiple software using a PairGrid.
    
    Parameters:
    - software_dfs: dict, keys are software names, values are DataFrames returned by extract_dpsi.

    Returns:
    - matplotlib.figure.Figure object.
    """
    plt.rcParams.update({
        'font.size': 11,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12
    })

    dfs_renamed = {
        name: df.rename(columns={'value': name})[['uniform_ID', name]]
        for name, df in software_dfs.items()
    }

    merged_df = None
    for df in dfs_renamed.values():
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on='uniform_ID', how='outer')

    data = merged_df.drop(columns='uniform_ID')

    def lower_scatter(x, y, **kws):
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() > 0:
            sns.scatterplot(x=x[mask], y=y[mask], alpha=0.6, s=15)

    def upper_corr(x, y, **kws):
        mask = np.isfinite(x) & np.isfinite(y)
        ax = plt.gca()
        if mask.sum() > 1:
            r, _ = pearsonr(x[mask], y[mask])
            ax.annotate(f"r = {r:.2f}\nn = {mask.sum()}",
                        xy=(.5, .5), xycoords='axes fraction',
                        ha='center', fontsize=11)
        else:
            ax.annotate("n < 2",
                        xy=(.5, .5), xycoords='axes fraction',
                        fontsize=11)

    def safe_hist(x, **kwargs):
        x = x[np.isfinite(x)]
        if len(x) > 0:
            plt.hist(x, **kwargs)
        else:
            plt.text(0.5, 0.5, "No valid data", ha='center', va='center', fontsize=12)
            plt.gca().set_xticks([])
            plt.gca().set_yticks([])

    g = sns.PairGrid(data, height=2)
    g.map_lower(lower_scatter)
    g.map_upper(upper_corr)
    g.map_diag(safe_hist)

    for ax in g.axes.flatten():
        if ax is not None:
            ax.title.set_fontsize(14)
            ax.xaxis.label.set_size(12)
            ax.yaxis.label.set_size(12)
            ax.tick_params(axis='both', labelsize=12)

    g.fig.set_size_inches(6.5, 6)
    g.fig.set_dpi(300)
    g.fig.tight_layout(rect=[0, 0, 1, 0.97])
    return g.fig


def plot_all_events_grid(event_dict, software_list, output_dir, sample_name, filename, ncols=2, figsize_scale=5):
    """
    Combine correlation plots of multiple events into a grid and save as an image.
    
    Parameters:
    - event_dict: dict, keys are event names, values are dicts {software_name: DataFrame}.
    - software_list: list of software names.
    - output_dir: str, directory to save the output image.
    - sample_name: str, sample identifier.
    - filename: str, output file name.
    - ncols: int, number of columns in the grid.
    - figsize_scale: int, scaling factor for figure size.
    """
    os.makedirs(os.path.join(output_dir, sample_name), exist_ok=True)

    n_events = len(event_dict)
    nrows = (n_events + ncols - 1) // ncols
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols,
                            figsize=(figsize_scale * ncols, figsize_scale * nrows))
    axs = axs.flatten()

    for i, (event_name, software_dfs) in enumerate(event_dict.items()):
        fig_i = plot_corr_dpsi(software_dfs)
        canvas = FigureCanvas(fig_i)
        canvas.draw()
        img = np.frombuffer(canvas.tostring_rgb(), dtype='uint8')
        img = img.reshape(fig_i.canvas.get_width_height()[::-1] + (3,))
        plt.close(fig_i)

        axs[i].imshow(img)
        axs[i].axis('off')
        axs[i].set_title(f"{event_name}", fontsize=18)

    for j in range(i + 1, len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    out_path = os.path.join(output_dir, sample_name, filename)
    plt.savefig(out_path, dpi=300)
    plt.close()
