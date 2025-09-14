# dse_plot_utils.py

import os
import math
import shutil
import warnings
from typing import Dict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from upsetplot import UpSet, from_contents
from PIL import Image


def plot_multiple_event_summaries(event_summary_dict: Dict[str, pd.DataFrame],
                                  output_dir: str,
                                  sample_name: str,
                                  filename: str):
    """
    Plot stacked bar charts for multiple DSE event summaries per software.
    
    event_summary_dict: {event_type: DataFrame with columns ['Software', 'up-regulate', 'down-regulate']}
    """
    num_events = len(event_summary_dict)
    if num_events == 0:
        return

    ncols = 2
    nrows = math.ceil(num_events / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6 * ncols, 5 * nrows), sharey=False)
    axes = axes.flatten() if num_events > 1 else [axes]

    for ax, (event, summary_df) in zip(axes, event_summary_dict.items()):
        df_plot = summary_df.set_index('Software')[['down-regulate', 'up-regulate']]
        plot = df_plot.plot(kind='bar', stacked=True, ax=ax, colormap='coolwarm', legend=False)

        for container in plot.containers:
            for bar in container:
                height = bar.get_height()
                if height > 0:
                    ax.annotate(f'{int(height)}',
                                xy=(bar.get_x() + bar.get_width() / 2, bar.get_y() + height / 2),
                                ha='center', va='center', fontsize=10, color='white')

        ax.set_title(f"{event}", fontsize=14)
        ax.set_xlabel("")
        ax.set_ylabel("DSE Event Count")

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, title='Regulation', loc='upper center', ncol=2)
    for ax in axes[num_events:]:
        ax.axis('off')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(os.path.join(output_dir, sample_name), exist_ok=True)
    save_path = os.path.join(output_dir, sample_name, filename)
    plt.savefig(save_path, dpi=300)
    plt.close()


def plot_venn_per_event_vennlib(event_dict: Dict[str, dict],
                                output_dir: str,
                                sample_name: str,
                                filename: str):
    """
    Plot Venn diagrams per DSE event showing overlap between software.
    """
    valid_events = {ev: data for ev, data in event_dict.items() if isinstance(data, dict) and len(data) > 1}
    num_events = len(valid_events)
    if num_events == 0:
        return

    ncols = 2
    nrows = math.ceil(num_events / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(8 * ncols, 6 * nrows))
    axes = axes.flatten() if num_events > 1 else [axes]

    for ax, (event, software_sets) in zip(axes, valid_events.items()):
        plt.sca(ax)
        try:
            venn(software_sets, ax=ax)
            ax.set_title(f'{event}', fontsize=14)
        except Exception as e:
            ax.text(0.5, 0.5, f"Error in {event}\n{e}", ha='center')
            ax.set_title(f'{event}')
            ax.axis('off')

    for ax in axes[num_events:]:
        ax.axis('off')

    plt.tight_layout()
    os.makedirs(os.path.join(output_dir, sample_name), exist_ok=True)
    save_path = os.path.join(output_dir, sample_name, filename)
    plt.savefig(save_path, dpi=300)
    plt.close()


def plot_upset_per_event_combined(event_dict: Dict[str, dict],
                                  output_dir: str,
                                  sample_name: str,
                                  filename: str = "combined_upset.png",
                                  ncols: int = 2):
    """
    Plot UpSet plots for multiple DSE events and combine them into a single image.
    """
    warnings.simplefilter(action='ignore', category=FutureWarning)
    temp_dir = os.path.join(output_dir, sample_name, "temp_upset_imgs")
    os.makedirs(temp_dir, exist_ok=True)
    image_paths = []

    for event, software_dict in event_dict.items():
        if len(software_dict) < 2:
            continue

        try:
            data = from_contents(software_dict)
            upset = UpSet(data, show_counts=True)

            fig = plt.figure(figsize=(6, 5))
            upset.plot(fig=fig)
            fig.suptitle(event, fontsize=14)

            save_path = os.path.join(temp_dir, f"{event}.png")
            fig.savefig(save_path, dpi=150, bbox_inches="tight")
            image_paths.append(save_path)
            plt.close(fig)
        except Exception as e:
            print(f"Error processing event {event}: {e}")
            continue

    if not image_paths:
        print("No valid events with ≥2 software for UpSet plotting.")
        return

    images = [Image.open(p) for p in image_paths]
    widths, heights = zip(*(img.size for img in images))
    max_width = max(widths)
    max_height = max(heights)

    nrows = (len(images) + ncols - 1) // ncols
    combined_image = Image.new("RGB", (max_width * ncols, max_height * nrows), "white")

    for idx, img in enumerate(images):
        row, col = divmod(idx, ncols)
        combined_image.paste(img, (col * max_width, row * max_height))

    os.makedirs(os.path.join(output_dir, sample_name), exist_ok=True)
    combined_path = os.path.join(output_dir, sample_name, filename)
    combined_image.save(combined_path)
    shutil.rmtree(temp_dir)
