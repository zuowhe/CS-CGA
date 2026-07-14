# -*- coding: utf-8 -*-
"""Plot averaged BIC convergence curves from per-run CSV files."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FuncFormatter


SCRIPT_DIR = Path(__file__).resolve().parent
PAPER_DIR = SCRIPT_DIR.parent
BIC_DATA_DIR = PAPER_DIR / "Result" / "BIC curve data"
BIC_OUTPUT_DIR = SCRIPT_DIR / "output" / "bic_convergence"

NETWORKS = ["Asia", "Insurance", "Water", "Alarm", "Hailfinder", "HeparII", "Win95pts", "Andes"]
DATASET_SIZES = [500, 1000, 3000]
ALGORITHMS = ["MAGA", "hybrid_SLA", "EKGA", "AESL_GA", "PSX", "MIGA", "CS-CGA"]

ITERATIONS = 200
X_VALUES = np.arange(ITERATIONS)
LINE_WIDTH = 2.5
MARKER_EVERY = 20
MARKER_SIZE = 8

CUSTOM_COLORS = {
    "MAGA": "#f08a5d",
    "hybrid_SLA": "#7ABBDB",
    "EKGA": "#84BA42",
    "AESL_GA": "#A5AEB7",
    "PSX": "#DBB428",
    "MIGA": "#DCA7EB",
    "CS-CGA": "#E34A33",
}

CUSTOM_MARKERS = {
    "MAGA": "p",
    "hybrid_SLA": "o",
    "EKGA": "s",
    "AESL_GA": "^",
    "PSX": "D",
    "MIGA": "v",
    "CS-CGA": "*",
}

DISPLAY_NAMES = {
    "MAGA": "MAGABN",
    "hybrid_SLA": "hybrid_SLA",
    "EKGA": "EKGA",
    "AESL_GA": "AESL_GA",
    "PSX": "PSX",
    "MIGA": "MIGA",
    "CS-CGA": "CS-CGA",
}

DATASET_LABELS = {
    ("Asia", 500): "AS-500",
    ("Asia", 1000): "AS-1000",
    ("Asia", 3000): "AS-3000",
    ("Insurance", 500): "INS-500",
    ("Insurance", 1000): "IN-1000",
    ("Insurance", 3000): "IN-3000",
    ("Water", 500): "WA-500",
    ("Water", 1000): "WA-1000",
    ("Water", 3000): "WA-3000",
    ("Alarm", 500): "AL-500",
    ("Alarm", 1000): "AL-1000",
    ("Alarm", 3000): "AL-3000",
    ("Hailfinder", 500): "HA-500",
    ("Hailfinder", 1000): "HA-1000",
    ("Hailfinder", 3000): "HA-3000",
    ("HeparII", 500): "HE-500",
    ("HeparII", 1000): "HE-1000",
    ("HeparII", 3000): "HE-3000",
    ("Win95pts", 500): "WI-500",
    ("Win95pts", 1000): "WI-1000",
    ("Win95pts", 3000): "WI-3000",
    ("Andes", 500): "AN-500",
    ("Andes", 1000): "AN-1000",
    ("Andes", 3000): "AN-3000",
}

DEFAULT_FONTS = {
    "base": 22,
    "title": 28,
    "axis_label": 26,
    "tick": 22,
    "legend": 20,
    "offset": 22,
}


def configure_plot_style(fonts=None):
    fonts = {**DEFAULT_FONTS, **(fonts or {})}
    plt.rcParams.update(
        {
            "font.size": fonts["base"],
            "axes.titlesize": fonts["title"],
            "axes.labelsize": fonts["axis_label"],
            "xtick.labelsize": fonts["tick"],
            "ytick.labelsize": fonts["tick"],
            "legend.fontsize": fonts["legend"],
        }
    )
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["axes.unicode_minus"] = False
    sns.set(style="whitegrid")
    plt.rcParams.update(
        {
            "axes.edgecolor": "#cccccc",
            "axes.linewidth": 1.0,
            "xtick.color": "black",
            "ytick.color": "black",
        }
    )
    return fonts


def parse_filename(filename):
    """Parse files like Asia1000_100_200_CS-CGA-01.csv."""
    if not filename.endswith(".csv"):
        return None

    name = filename[:-4]

    try:
        # Algorithm names may contain hyphens, so split only the trailing run index.
        main_part, index_part = name.rsplit("-", 1)
        if len(index_part) != 2 or not index_part.isdigit():
            return None

        parts = main_part.split("_")
        if len(parts) < 4:
            return None

        net_data_part = parts[0]
        pop_size_str = parts[1]
        iterations_str = parts[2]
        algorithm = "_".join(parts[3:])

        data_size_str = ""
        i = len(net_data_part) - 1
        while i >= 0 and net_data_part[i].isdigit():
            data_size_str = net_data_part[i] + data_size_str
            i -= 1

        if not data_size_str:
            return None

        net_name = net_data_part[: i + 1]
        if net_name not in NETWORKS:
            return None

        return {
            "net_name": net_name,
            "data_size": int(data_size_str),
            "pop_size": int(pop_size_str),
            "iterations": int(iterations_str),
            "algorithm": algorithm,
            "index": index_part,
        }
    except (ValueError, IndexError):
        return None


def diagnose_data_availability():
    """Print a summary of available CSV files before plotting."""
    print("Checking BIC curve data...")
    print(f"Data directory: {BIC_DATA_DIR}")

    if not BIC_DATA_DIR.exists():
        print(f"Data directory does not exist: {BIC_DATA_DIR}")
        return

    network_counts = {net: 0 for net in NETWORKS}

    for algo in ALGORITHMS:
        algo_path = BIC_DATA_DIR / algo
        if not algo_path.exists():
            print(f"Missing algorithm directory: {algo_path}")
            continue

        print(f"\nAlgorithm: {algo}")
        files_found = 0
        unmatched_files = []
        network_files = {net: [] for net in NETWORKS}

        for file_path in algo_path.iterdir():
            if file_path.suffix.lower() != ".csv":
                continue

            parsed = parse_filename(file_path.name)
            if parsed:
                files_found += 1
                net_name = parsed["net_name"]
                network_counts[net_name] += 1
                network_files[net_name].append((file_path.name, parsed["data_size"], parsed["iterations"]))
            else:
                unmatched_files.append(file_path.name)

        print(f"  Parsed CSV files: {files_found}")
        for net in NETWORKS:
            if network_files[net]:
                print(f"  {net}: {len(network_files[net])} files")
                for filename, size, iterations in network_files[net][:3]:
                    print(f"    - {filename} (size={size}, iterations={iterations})")

        if unmatched_files:
            print(f"  Unmatched files: {len(unmatched_files)}")
            for filename in unmatched_files[:5]:
                print(f"    - {filename}")

    print("\nOverall parsed file counts:")
    for net in NETWORKS:
        print(f"  {net}: {network_counts[net]}")


def read_bic_curve(file_path):
    df = pd.read_csv(file_path, header=None)
    if df.empty or df.shape[1] == 0:
        return None

    bic_values = df.iloc[:, 0].to_numpy(dtype=float)
    if len(bic_values) == 0:
        return None

    if len(bic_values) < ITERATIONS:
        return np.concatenate([bic_values, np.full(ITERATIONS - len(bic_values), bic_values[-1])])

    return bic_values[:ITERATIONS]


def use_scientific_y_axis(ax, fonts, tick_count=5, inset_ratio=0.08):
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    inner_min = y_min + y_range * inset_ratio
    inner_max = y_max - y_range * inset_ratio
    tick_values = np.linspace(inner_min, inner_max, tick_count)
    ax.set_yticks(tick_values)

    max_abs_tick = max(abs(value) for value in tick_values)
    exponent = int(np.floor(np.log10(max_abs_tick))) if max_abs_tick > 0 else 0
    scale = 10**exponent if exponent else 1

    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / scale:.2f}"))
    ax.yaxis.get_offset_text().set_visible(False)
    if exponent:
        ax.text(
            0.02,
            0.98,
            rf"$\times 10^{{{exponent}}}$",
            transform=ax.transAxes,
            fontsize=fonts["offset"],
            horizontalalignment="left",
            verticalalignment="top",
        )


def process_and_plot(fonts):
    """Average each algorithm's BIC curves and save convergence PDFs."""
    output_dir = BIC_OUTPUT_DIR
    output_dir.mkdir(parents=True, exist_ok=True)
    fallback_palette = sns.color_palette("husl", len(ALGORITHMS))

    for network in NETWORKS:
        for size in DATASET_SIZES:
            fig, ax = plt.subplots(figsize=(8.5, 6.4))
            for spine in ax.spines.values():
                spine.set_color("#cccccc")
                spine.set_linewidth(1.0)
            title = f"{network}{size}"
            plotted = False
            debug_info = {algo: 0 for algo in ALGORITHMS}

            print(f"\nProcessing: {title}")

            for algo in ALGORITHMS:
                folder_path = BIC_DATA_DIR / algo
                if not folder_path.exists():
                    print(f"  Missing directory for {algo}: {folder_path}")
                    continue

                bic_curves = []

                for file_path in folder_path.iterdir():
                    if file_path.suffix.lower() != ".csv":
                        continue

                    parsed = parse_filename(file_path.name)
                    if not parsed:
                        continue

                    if (
                        parsed["net_name"] != network
                        or parsed["data_size"] != size
                        or parsed["algorithm"] != algo
                        or parsed["iterations"] != ITERATIONS
                    ):
                        continue

                    try:
                        bic_curve = read_bic_curve(file_path)
                    except Exception as exc:
                        print(f"  Failed to read {file_path.name}: {exc}")
                        continue

                    if bic_curve is not None:
                        bic_curves.append(bic_curve)

                debug_info[algo] = len(bic_curves)

                if len(bic_curves) >= 10:
                    mean_curve = np.mean(np.array(bic_curves), axis=0)
                    color = CUSTOM_COLORS.get(algo, fallback_palette[ALGORITHMS.index(algo)])
                    marker = CUSTOM_MARKERS.get(algo, "o")
                    label = DISPLAY_NAMES.get(algo, algo)

                    ax.plot(
                        X_VALUES,
                        mean_curve,
                        label=label,
                        color=color,
                        linewidth=LINE_WIDTH,
                        marker=marker,
                        markevery=MARKER_EVERY,
                        markersize=MARKER_SIZE,
                        markerfacecolor=color,
                    )
                    plotted = True
                elif bic_curves:
                    print(f"  {algo}: found {len(bic_curves)} files, fewer than the required 10")
                else:
                    print(f"  {algo}: no matching files")

            print(f"  {title} file counts: {debug_info}")

            if plotted:
                dataset_label = DATASET_LABELS[(network, size)]
                ax.set_title(f"BIC Convergence - {dataset_label}", fontsize=fonts["title"], pad=20)
                ax.set_xlabel("Iteration", fontsize=fonts["axis_label"])
                ax.set_ylabel("Mean BIC Score", fontsize=fonts["axis_label"])
                ax.tick_params(axis="both", labelsize=fonts["tick"])
                use_scientific_y_axis(ax, fonts)
                ax.legend(loc="lower right", fontsize=fonts["legend"], frameon=True, fancybox=True, shadow=False)
                ax.grid(False)
                ax.set_xlim(0, ITERATIONS - 1)
                plt.margins(x=0)
                plt.tight_layout()

                output_path = output_dir / f"{title}_convergence.pdf"
                plt.savefig(output_path, dpi=300, bbox_inches="tight")
                print(f"  Saved: {output_path}")
            else:
                print(f"  Skipped: {title} has no sufficient data")

            plt.close(fig)


def main(fonts=None):
    fonts = configure_plot_style(fonts)
    diagnose_data_availability()
    process_and_plot(fonts)
    print("Done.")


if __name__ == "__main__":
    main()
