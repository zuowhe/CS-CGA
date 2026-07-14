# -*- coding: utf-8 -*-
"""Plot F1 curves for adaptive threshold embedding in hybrid SLA."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_FILE = SCRIPT_DIR / "data" / "adaptive_threshold_sla.csv"
OUTPUT_DIR = SCRIPT_DIR / "output" / "adaptive_threshold_sla"

DEFAULT_FONTS = {
    "base": 20,
    "title": 26,
    "axis_label": 26,
    "tick": 23,
    "legend": 19,
}


def configure_plot_style(fonts=None):
    fonts = {**DEFAULT_FONTS, **(fonts or {})}
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams.update(
        {
            "font.size": fonts["base"],
            "axes.titlesize": fonts["title"],
            "axes.labelsize": fonts["axis_label"],
            "xtick.labelsize": fonts["tick"],
            "ytick.labelsize": fonts["tick"],
            "legend.fontsize": fonts["legend"],
            "axes.edgecolor": "black",
            "axes.linewidth": 1.0,
            "xtick.color": "black",
            "ytick.color": "black",
        }
    )
    return fonts


def load_threshold_data(filename):
    adaptive_f1 = {}
    fixed_points = defaultdict(list)

    with filename.open("r", encoding="utf-8", newline="") as file_obj:
        reader = csv.DictReader(file_obj)
        for row in reader:
            dataset = row["dataset"]
            series = row["series"]
            mean = float(row["mean"])

            if series == "adaptive":
                adaptive_f1[dataset] = mean
            elif series == "fixed":
                fixed_points[dataset].append((float(row["alpha"]), mean))

    for dataset in fixed_points:
        fixed_points[dataset].sort(key=lambda item: item[0])

    return adaptive_f1, fixed_points


def add_top_legend_space(ax, values):
    y_min = min(values)
    y_max = max(values)
    y_range = y_max - y_min
    if y_range <= 0:
        y_range = max(abs(y_max) * 0.03, 0.02)

    lower_padding = y_range * 0.12
    upper_padding = y_range * 0.28
    ax.set_ylim(y_min - lower_padding, y_max + upper_padding)


def format_y_axis(ax, tick_count=6, inset_ratio=0.08):
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    inner_min = y_min + y_range * inset_ratio
    inner_max = y_max - y_range * inset_ratio
    step = (inner_max - inner_min) / (tick_count - 1)
    ax.set_yticks([inner_min + step * index for index in range(tick_count)])
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))


def plot_dataset(dataset, adaptive_value, fixed_points, output_dir, fonts):
    alpha_values = [point[0] for point in fixed_points]
    f1_values = [point[1] for point in fixed_points]

    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    for spine in ax.spines.values():
        spine.set_color("black")
        spine.set_linewidth(1.0)
    ax.plot(alpha_values, f1_values, marker="o", label=r"Fixed $\alpha$")
    ax.fill_between(alpha_values, adaptive_value, f1_values, color="skyblue", alpha=0.4)
    ax.axhline(
        y=adaptive_value,
        color="r",
        linestyle="--",
        linewidth=2,
        label=r"Adaptive $\alpha$",
    )

    ax.set_xlabel(r"The value of $\alpha$", fontsize=fonts["axis_label"])
    ax.set_ylabel("F1 Score", fontsize=fonts["axis_label"])
    ax.set_xscale("log")
    ax.tick_params(axis="both", labelsize=fonts["tick"])
    add_top_legend_space(ax, f1_values + [adaptive_value])
    format_y_axis(ax)
    ax.legend(
        fontsize=fonts["legend"],
        loc="upper center",
        ncol=2,
        frameon=True,
    )
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.tight_layout()

    output_path = output_dir / f"{dataset}_SLA_F1.pdf"
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output_path}")


def main(fonts=None):
    fonts = configure_plot_style(fonts)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    adaptive_f1, fixed_points = load_threshold_data(DATA_FILE)

    for dataset, points in fixed_points.items():
        plot_dataset(dataset, adaptive_f1[dataset], points, OUTPUT_DIR, fonts)

    print("Done.")


if __name__ == "__main__":
    main()
