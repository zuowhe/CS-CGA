# -*- coding: utf-8 -*-
"""Plot running-time ablation bars from an external CSV file."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_FILE = SCRIPT_DIR / "data" / "time_ablation_runtime.csv"
OUTPUT_DIR = SCRIPT_DIR / "output" / "time_ablation"

STRATEGIES = ["No-Cache", "BNT-Cache", "GB-Cache", "HSM"]
COLORS = ["#8AD1F2", "#6DE87D", "#DDBFDE", "#E54046"]

DEFAULT_FONTS = {
    "base": 22,
    "axes_label": 29,
    "x_tick": 30,
    "y_tick": 32,
    "legend": 29,
    "y_label": 34,
}


def configure_plot_style(fonts=None):
    fonts = {**DEFAULT_FONTS, **(fonts or {})}
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams.update(
        {
            "font.size": fonts["base"],
            "axes.labelsize": fonts["axes_label"],
            "xtick.labelsize": fonts["x_tick"],
            "ytick.labelsize": fonts["y_tick"],
            "legend.fontsize": fonts["legend"],
            "axes.edgecolor": "black",
            "axes.linewidth": 1.0,
            "xtick.color": "black",
            "ytick.color": "black",
        }
    )
    return fonts


def load_runtime_data(filename):
    rows = []
    with filename.open("r", encoding="utf-8", newline="") as file_obj:
        reader = csv.DictReader(file_obj)
        for row in reader:
            rows.append(
                {
                    "dataset": row["dataset"],
                    "strategy": row["strategy"],
                    "mean": float(row["mean"]),
                    "std": float(row["std"]),
                }
            )
    return rows


def split_rows(rows):
    datasets = []
    for row in rows:
        if row["dataset"] not in datasets:
            datasets.append(row["dataset"])

    split_index = datasets.index("HA-500")
    before = set(datasets[:split_index])
    after = set(datasets[split_index:])

    rows_before = [row for row in rows if row["dataset"] in before]
    rows_after = [row for row in rows if row["dataset"] in after]
    return rows_before, rows_after


def matrix_from_rows(rows):
    datasets = []
    for row in rows:
        if row["dataset"] not in datasets:
            datasets.append(row["dataset"])

    means = np.zeros((len(datasets), len(STRATEGIES)))
    stds = np.zeros((len(datasets), len(STRATEGIES)))
    dataset_index = {dataset: idx for idx, dataset in enumerate(datasets)}
    strategy_index = {strategy: idx for idx, strategy in enumerate(STRATEGIES)}

    for row in rows:
        i = dataset_index[row["dataset"]]
        j = strategy_index[row["strategy"]]
        means[i, j] = row["mean"]
        stds[i, j] = row["std"]

    return datasets, means, stds


def plot_and_save(rows, filename, fonts):
    datasets, means, stds = matrix_from_rows(rows)
    fig, ax = plt.subplots(figsize=(16, 10))
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.grid(False)
    ax.tick_params(axis="y", which="both", length=5, direction="out", labelsize=fonts["y_tick"])

    n_datasets = len(datasets)
    n_strategies = len(STRATEGIES)
    bar_width = 0.8 / n_strategies
    index = np.arange(n_datasets)

    for strategy_idx, strategy in enumerate(STRATEGIES):
        ax.bar(
            index + strategy_idx * bar_width,
            means[:, strategy_idx],
            bar_width,
            yerr=stds[:, strategy_idx],
            capsize=3,
            label=strategy,
            color=COLORS[strategy_idx],
            alpha=0.8,
            zorder=2,
        )

    ax.set_ylabel("Running Time (seconds)", fontsize=fonts["y_label"])
    ax.set_xticks(index + bar_width * (n_strategies - 1) / 2)
    ax.set_xticklabels(datasets, rotation=45, ha="right", fontsize=fonts["x_tick"])
    ax.legend(fontsize=fonts["legend"], loc="upper left")
    plt.tight_layout()

    output_path = OUTPUT_DIR / filename
    with PdfPages(output_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main(fonts=None):
    fonts = configure_plot_style(fonts)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = load_runtime_data(DATA_FILE)
    rows_before, rows_after = split_rows(rows)
    plot_and_save(rows_before, "running_time_before_HA500.pdf", fonts)
    plot_and_save(rows_after, "running_time_from_HA500.pdf", fonts)
    print("Done.")


if __name__ == "__main__":
    main()
