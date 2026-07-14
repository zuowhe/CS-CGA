# -*- coding: utf-8 -*-
"""Plot normalized F1 and SHD radar charts from an external CSV file."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "data"
OUTPUT_DIR = SCRIPT_DIR / "output" / "radar_charts"
METRIC_DATA_FILE = DATA_DIR / "radar_metrics.csv"

DATASETS = ["Asia", "Insurance", "Water", "Alarm", "Hailfinder", "HeparII", "Win95pts", "AND"]
SAMPLE_SIZES = [500, 1000, 3000]
METHODS = ["SR-Fix", "SR-PCR", "SM-Fix", "SM-PCR", "DR-PCR", "CS-CGA"]
COLORS = ["#a2738c", "#6eb6ff", "#66c6ba", "#ffd460", "#3A5FCD", "#ea5455"]

DEFAULT_FONTS = {
    "base": 22,
    "tick": 22,
    "legend": 22,
}


def configure_plot_style(fonts=None):
    fonts = {**DEFAULT_FONTS, **(fonts or {})}
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams.update(
        {
            "font.size": fonts["base"],
            "xtick.labelsize": fonts["tick"],
            "ytick.labelsize": fonts["tick"],
            "legend.fontsize": fonts["legend"],
            "axes.edgecolor": "#cccccc",
            "axes.linewidth": 1.0,
            "xtick.color": "black",
            "ytick.color": "black",
        }
    )
    return fonts


def load_metric_data(filename):
    rows = {}
    with filename.open("r", encoding="utf-8", newline="") as file_obj:
        reader = csv.DictReader(file_obj)
        for row in reader:
            metric = row["metric"]
            dataset = row["dataset"]
            sample_size = int(row["sample_size"])
            values = [float(row[method]) for method in METHODS]
            rows[(metric, dataset, sample_size)] = values
    return rows


def metric_matrix(rows, metric, sample_size):
    values = []
    for dataset in DATASETS:
        values.append(rows[(metric, dataset, sample_size)])
    return np.array(values, dtype=float)


def normalize_scores(scores, higher_is_better=True):
    normalized = np.zeros_like(scores)
    for row_idx in range(scores.shape[0]):
        row = scores[row_idx]
        min_value = row.min()
        max_value = row.max()

        if max_value == min_value:
            normalized[row_idx, :] = 1.0
        elif higher_is_better:
            normalized[row_idx, :] = (row - min_value) / (max_value - min_value)
        else:
            normalized[row_idx, :] = (max_value - row) / (max_value - min_value)

    return normalized


def build_plot_configs(rows):
    configs = []
    for sample_size in SAMPLE_SIZES:
        f1_scores = metric_matrix(rows, "F1", sample_size)
        configs.append((f"F1_{sample_size}", normalize_scores(f1_scores, higher_is_better=True)))

    for sample_size in SAMPLE_SIZES:
        shd_scores = metric_matrix(rows, "SHD", sample_size)
        configs.append((f"SHD_{sample_size}", normalize_scores(shd_scores, higher_is_better=False)))

    return configs


def plot_radar_chart(title, data, output_dir, fonts):
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={"projection": "polar"})
    for spine in ax.spines.values():
        spine.set_color("#cccccc")
        spine.set_linewidth(1.0)

    angles = np.linspace(0, 2 * np.pi, len(DATASETS), endpoint=False).tolist()
    angles += angles[:1]

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetagrids(np.degrees(angles[:-1]), DATASETS)
    ax.set_ylim(-0.2, 1)
    ax.tick_params(axis="y", labelsize=fonts["tick"], labelcolor="black")
    ax.tick_params(axis="x", labelsize=fonts["tick"], labelcolor="black")

    for method_idx, method in enumerate(METHODS):
        values = data[:, method_idx].tolist()
        values += values[:1]
        color = COLORS[method_idx]
        alpha_fill = 0.25 if method_idx == len(METHODS) - 1 else 0.10
        ax.plot(angles, values, linewidth=2, label=method, color=color)
        ax.fill(angles, values, alpha=alpha_fill, color=color)

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=3,
        prop={"family": "Times New Roman", "size": fonts["legend"]},
    )
    plt.tight_layout()

    output_path = output_dir / f"{title}_radar.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main(fonts=None):
    fonts = configure_plot_style(fonts)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_metric_data(METRIC_DATA_FILE)
    for title, data in build_plot_configs(rows):
        plot_radar_chart(title, data, OUTPUT_DIR, fonts)

    print("Done.")


if __name__ == "__main__":
    main()
