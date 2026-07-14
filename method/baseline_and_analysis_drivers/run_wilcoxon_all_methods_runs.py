from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


DATASETS = ["Asia", "INS", "Water", "Alarm", "Hailfinder", "HEPAR", "Win95pts", "AND"]
SAMPLE_SIZES = [500, 1000, 3000]
REFERENCE_ALGORITHM = "CS-CGA"
METRICS = {
    "F1": {"column": 1, "higher_is_better": True},
    "SHD": {"column": 5, "higher_is_better": False},
    "Runtime": {"column": 7, "higher_is_better": False},
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Run trial-level Wilcoxon tests for all current manuscript baselines."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "my paper" / "Result" / "Wilcoxon_all_runs",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def comparison_sources(root: Path) -> dict[str, Path]:
    paper_results = root / "my paper" / "Result"
    return {
        "PSX": paper_results / "PSX",
        "EKGA-BN": paper_results / "EKGA",
        "AESL-GA": paper_results / "AESL",
        "hybrid-SLA": paper_results / "SLA",
        "MIGA": paper_results / "MIGA",
        "MAGABN": paper_results / "MAGABN",
        "GTT": root / "results" / "GTT",
        "BS": root / "results" / "BS",
        "IGES": root / "results" / "IGES",
        "saiyanH": root / "results" / "saiyanH",
    }


def parse_instance(path: Path) -> tuple[str, int] | None:
    for dataset in sorted(DATASETS, key=len, reverse=True):
        if not path.name.startswith(dataset):
            continue
        match = re.match(rf"{re.escape(dataset)}(\d+)_100_200_", path.name)
        if match is None:
            return None
        sample_size = int(match.group(1))
        if sample_size in SAMPLE_SIZES:
            return dataset, sample_size
    return None


def build_index(directory: Path) -> dict[tuple[str, int], Path]:
    index: dict[tuple[str, int], Path] = {}
    for path in sorted(directory.glob("*.csv")):
        instance = parse_instance(path)
        if instance is not None:
            index.setdefault(instance, path)
    return index


def try_float(value: str) -> float | None:
    try:
        number = float(value.strip())
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def read_trials(path: Path) -> dict[int, dict[str, float]]:
    """Read trials 1--10, retaining the final record for any repeated trial ID."""
    trials: dict[int, dict[str, float]] = {}
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        for raw_line in handle:
            cells = [cell.strip() for cell in raw_line.strip().split(",")]
            if len(cells) <= max(spec["column"] for spec in METRICS.values()):
                continue
            trial_value = try_float(cells[0])
            if trial_value is None or not trial_value.is_integer():
                continue
            trial = int(trial_value)
            if trial not in range(1, 11):
                continue
            values = {
                metric: try_float(cells[spec["column"]])
                for metric, spec in METRICS.items()
            }
            if all(value is not None for value in values.values()):
                trials[trial] = {metric: float(value) for metric, value in values.items()}
    return trials


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    position = 0
    while position < len(values):
        end = position + 1
        while end < len(values) and values[order[end]] == values[order[position]]:
            end += 1
        rank = (position + 1 + end) / 2.0
        for rank_index in range(position, end):
            ranks[order[rank_index]] = rank
        position = end
    return ranks


def normal_survival(value: float) -> float:
    """Numerically stable upper-tail probability for the standard normal law."""
    return 0.5 * math.erfc(value / math.sqrt(2.0))


def wilcoxon_signed_rank(improvements: list[float]) -> dict[str, float | int | str]:
    nonzero = [value for value in improvements if abs(value) > 1e-12]
    if not nonzero:
        return {"n_nonzero": 0, "w_plus": 0.0, "w_minus": 0.0, "p_raw": 1.0}

    ranks = average_ranks([abs(value) for value in nonzero])
    w_plus = sum(rank for rank, value in zip(ranks, nonzero) if value > 0)
    w_minus = sum(rank for rank, value in zip(ranks, nonzero) if value < 0)
    mean_rank_sum = sum(ranks) / 2.0
    variance = sum(rank * rank for rank in ranks) / 4.0
    if variance == 0:
        p_raw = 1.0
    else:
        correction = 0.5 if w_plus > mean_rank_sum else -0.5 if w_plus < mean_rank_sum else 0.0
        z_value = (w_plus - mean_rank_sum - correction) / math.sqrt(variance)
        p_raw = min(1.0, 2.0 * normal_survival(abs(z_value)))
    return {
        "n_nonzero": len(nonzero),
        "w_plus": w_plus,
        "w_minus": w_minus,
        "p_raw": p_raw,
    }


def holm_adjust(rows: list[dict[str, object]]) -> None:
    ordered = sorted(rows, key=lambda row: float(row["p_raw"]))
    total_tests = len(ordered)
    running_max = 0.0
    for index, row in enumerate(ordered):
        adjusted = min(1.0, (total_tests - index) * float(row["p_raw"]))
        running_max = max(running_max, adjusted)
        row["p_holm"] = running_max


def significance_label(p_value: float, alpha: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < alpha:
        return "*"
    return "ns"


def format_p(value: float) -> str:
    return f"{value:.2e}" if value < 1e-4 else f"{value:.4f}"


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_latex_table(path: Path, rows: list[dict[str, object]]) -> None:
    by_method = {
        method: {row["metric"]: row for row in rows if row["method"] == method}
        for method in comparison_sources(Path(__file__).resolve().parents[1])
    }
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\\begin{table}[htbp]\n")
        handle.write("\\centering\n")
        handle.write(
            "\\caption{Wilcoxon signed-rank tests based on matched ten-run results for CS-CGA against all compared methods. Holm correction is applied across all 30 method-metric tests.}\n"
        )
        handle.write("\\label{tab:wilcoxon_test}\n")
        handle.write("\\small\n")
        handle.write("\\begin{tabular}{lcccc}\n")
        handle.write("\\toprule\n")
        handle.write("Compared method & Pairs & F1 Holm-$p$ & SHD Holm-$p$ & Runtime Holm-$p$ \\\\" + "\n")
        handle.write("\\midrule\n")
        for method, metric_rows in by_method.items():
            pairs = next(iter(metric_rows.values()))["n_pairs"]
            cells = []
            for metric in METRICS:
                row = metric_rows[metric]
                sign = "+" if row["wins"] > row["losses"] else "-"
                cells.append(f"{format_p(float(row['p_holm']))} ({row['significance']}, {sign})")
            handle.write(f"{method} & {pairs} & " + " & ".join(cells) + " \\\\" + "\n")
        handle.write("\\bottomrule\n")
        handle.write("\\end{tabular}\n")
        handle.write("\\\\\n")
        handle.write(
            "\\footnotesize{Each pair comprises the same trial index under one dataset-sample setting. Positive signs indicate better CS-CGA performance, i.e., higher F1 or lower SHD/runtime.}\n"
        )
        handle.write("\\end{table}\n")


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_index = build_index(root / "my paper" / "Result" / REFERENCE_ALGORITHM)
    source_indices = {method: build_index(path) for method, path in comparison_sources(root).items()}
    summary_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []

    ordered_instances = sorted(
        reference_index,
        key=lambda instance: (DATASETS.index(instance[0]), SAMPLE_SIZES.index(instance[1])),
    )
    for method, comparison_index in source_indices.items():
        improvements = {metric: [] for metric in METRICS}
        wins = {metric: 0 for metric in METRICS}
        ties = {metric: 0 for metric in METRICS}
        losses = {metric: 0 for metric in METRICS}
        settings = 0
        for dataset, sample_size in ordered_instances:
            reference_path = reference_index[(dataset, sample_size)]
            comparison_path = comparison_index.get((dataset, sample_size))
            if comparison_path is None:
                missing_rows.append(
                    {"method": method, "dataset": dataset, "sample_size": sample_size, "reason": "missing_file"}
                )
                continue
            reference_trials = read_trials(reference_path)
            comparison_trials = read_trials(comparison_path)
            trial_ids = sorted(set(reference_trials) & set(comparison_trials))
            if trial_ids != list(range(1, 11)):
                missing_rows.append(
                    {
                        "method": method,
                        "dataset": dataset,
                        "sample_size": sample_size,
                        "reason": f"matched_trials_{len(trial_ids)}",
                    }
                )
                continue
            settings += 1
            for trial in trial_ids:
                for metric, spec in METRICS.items():
                    reference_value = reference_trials[trial][metric]
                    comparison_value = comparison_trials[trial][metric]
                    improvement = (
                        reference_value - comparison_value
                        if spec["higher_is_better"]
                        else comparison_value - reference_value
                    )
                    improvements[metric].append(improvement)
                    if improvement > 1e-12:
                        wins[metric] += 1
                    elif improvement < -1e-12:
                        losses[metric] += 1
                    else:
                        ties[metric] += 1
                    detail_rows.append(
                        {
                            "method": method,
                            "dataset": dataset,
                            "sample_size": sample_size,
                            "trial": trial,
                            "metric": metric,
                            "cscga_value": reference_value,
                            "comparison_value": comparison_value,
                            "cscga_improvement": improvement,
                        }
                    )
        for metric, values in improvements.items():
            test = wilcoxon_signed_rank(values)
            summary_rows.append(
                {
                    "method": method,
                    "metric": metric,
                    "matched_settings": settings,
                    "n_pairs": len(values),
                    "n_nonzero": test["n_nonzero"],
                    "wins": wins[metric],
                    "ties": ties[metric],
                    "losses": losses[metric],
                    "w_plus": test["w_plus"],
                    "w_minus": test["w_minus"],
                    "p_raw": test["p_raw"],
                }
            )

    holm_adjust(summary_rows)
    for row in summary_rows:
        row["significance"] = significance_label(float(row["p_holm"]), args.alpha)

    summary_fields = [
        "method", "metric", "matched_settings", "n_pairs", "n_nonzero", "wins", "ties", "losses",
        "w_plus", "w_minus", "p_raw", "p_holm", "significance",
    ]
    detail_fields = [
        "method", "dataset", "sample_size", "trial", "metric", "cscga_value", "comparison_value", "cscga_improvement",
    ]
    missing_fields = ["method", "dataset", "sample_size", "reason"]
    write_csv(output_dir / "wilcoxon_all_methods_trial_level.csv", detail_rows, detail_fields)
    write_csv(output_dir / "wilcoxon_all_methods_summary.csv", summary_rows, summary_fields)
    write_csv(output_dir / "wilcoxon_all_methods_missing.csv", missing_rows, missing_fields)
    write_latex_table(output_dir / "wilcoxon_all_methods_table.tex", summary_rows)
    print(f"Wrote: {output_dir / 'wilcoxon_all_methods_summary.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_all_methods_trial_level.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_all_methods_missing.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_all_methods_table.tex'}")


if __name__ == "__main__":
    main()
