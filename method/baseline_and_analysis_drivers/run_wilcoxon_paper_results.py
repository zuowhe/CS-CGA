from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


DATASETS = ["Asia", "INS", "Water", "Alarm", "Hailfinder", "HEPAR", "Win95pts", "AND"]
SAMPLE_SIZES = [500, 1000, 3000]
REFERENCE_ALGO = "CS-CGA"
COMPARISON_ALGOS = {
    "PSX": "PSX",
    "EKGA-BN": "EKGA",
    "AESL-GA": "AESL",
    "hybrid-SLA": "SLA",
    "MIGA": "MIGA",
    "MAGABN": "MAGABN",
}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    default_result_dir = root / "my paper" / "Result"
    default_output_dir = default_result_dir / "Wilcoxon"
    parser = argparse.ArgumentParser(
        description="Run Wilcoxon signed-rank tests for paper result CSV files."
    )
    parser.add_argument("--result-dir", type=Path, default=default_result_dir)
    parser.add_argument("--output-dir", type=Path, default=default_output_dir)
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Significance level after Holm correction.",
    )
    parser.add_argument(
        "--unit",
        choices=["instance-mean", "run"],
        default="instance-mean",
        help=(
            "Statistical unit for the Wilcoxon test. Use instance-mean for one "
            "paired value per dataset/sample-size setting, or run for all matched runs."
        ),
    )
    return parser.parse_args()


def parse_instance(path: Path) -> tuple[str, int] | None:
    name = path.name
    for dataset in sorted(DATASETS, key=len, reverse=True):
        if not name.startswith(dataset):
            continue
        rest = name[len(dataset) :]
        match = re.match(r"(\d+)_100_200_", rest)
        if not match:
            return None
        sample_size = int(match.group(1))
        if sample_size not in SAMPLE_SIZES:
            return None
        return dataset, sample_size
    return None


def build_file_index(algo_dir: Path) -> dict[tuple[str, int], Path]:
    index: dict[tuple[str, int], Path] = {}
    if not algo_dir.exists():
        return index

    for path in sorted(algo_dir.glob("*.csv")):
        lower_name = path.name.lower()
        if "summary" in lower_name or "merged" in lower_name:
            continue
        instance = parse_instance(path)
        if instance is None:
            continue
        index.setdefault(instance, path)
    return index


def try_float(value: str) -> float | None:
    try:
        return float(value.strip())
    except ValueError:
        return None


def read_f1_shd(path: Path) -> dict[str, list[float]]:
    f1_values: list[float] = []
    shd_values: list[float] = []
    with path.open("r", encoding="utf-8-sig", errors="ignore", newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            cells = [cell.strip() for cell in line.split(",")]
            if len(cells) < 6:
                continue
            iteration = try_float(cells[0])
            if iteration is None:
                continue
            f1 = try_float(cells[1])
            shd = try_float(cells[5])
            if f1 is None or shd is None:
                continue
            f1_values.append(f1)
            shd_values.append(shd)
    return {"F1": f1_values, "SHD": shd_values}


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    ranks = [0.0] * len(values)
    pos = 0
    while pos < len(values):
        end = pos + 1
        current = values[order[pos]]
        while end < len(values) and values[order[end]] == current:
            end += 1
        avg_rank = (pos + 1 + end) / 2.0
        for k in range(pos, end):
            ranks[order[k]] = avg_rank
        pos = end
    return ranks


def exact_signed_rank_p_value(scaled_ranks: list[int], observed_w_plus: int) -> float:
    total_rank_sum = sum(scaled_ranks)
    counts = [0] * (total_rank_sum + 1)
    counts[0] = 1
    max_seen = 0
    for rank in scaled_ranks:
        next_counts = counts[:]
        for score in range(max_seen + 1):
            count = counts[score]
            if count:
                next_counts[score + rank] += count
        counts = next_counts
        max_seen += rank

    low_tail = min(observed_w_plus, total_rank_sum - observed_w_plus)
    high_tail = total_rank_sum - low_tail
    tail_count = sum(counts[: low_tail + 1]) + sum(counts[high_tail:])
    total_count = 2 ** len(scaled_ranks)
    return min(1.0, tail_count / total_count)


def normal_cdf(value: float) -> float:
    return 0.5 * (1.0 + math.erf(value / math.sqrt(2.0)))


def normal_signed_rank_p_value(ranks: list[float], w_plus: float) -> float:
    mean = sum(ranks) / 2.0
    variance = sum(rank * rank for rank in ranks) / 4.0
    if variance <= 0:
        return 1.0
    correction = 0.5 if w_plus > mean else -0.5 if w_plus < mean else 0.0
    z = (w_plus - mean - correction) / math.sqrt(variance)
    return min(1.0, 2.0 * (1.0 - normal_cdf(abs(z))))


def wilcoxon_signed_rank(improvements: list[float]) -> dict[str, float | int | str]:
    nonzero = [value for value in improvements if abs(value) > 1e-12]
    n = len(nonzero)
    if n == 0:
        return {
            "n_nonzero": 0,
            "w_plus": 0.0,
            "w_minus": 0.0,
            "p_raw": 1.0,
            "p_method": "all_zero",
        }

    abs_values = [abs(value) for value in nonzero]
    ranks = average_ranks(abs_values)
    w_plus = sum(rank for rank, value in zip(ranks, nonzero) if value > 0)
    w_minus = sum(rank for rank, value in zip(ranks, nonzero) if value < 0)
    scaled_ranks = [int(round(rank * 2)) for rank in ranks]
    observed_scaled_w_plus = int(round(w_plus * 2))

    if sum(scaled_ranks) <= 200000:
        p_raw = exact_signed_rank_p_value(scaled_ranks, observed_scaled_w_plus)
        p_method = "exact_signed_rank_distribution"
    else:
        p_raw = normal_signed_rank_p_value(ranks, w_plus)
        p_method = "normal_approximation"

    return {
        "n_nonzero": n,
        "w_plus": w_plus,
        "w_minus": w_minus,
        "p_raw": p_raw,
        "p_method": p_method,
    }


def holm_adjust(rows: list[dict[str, object]], metric: str) -> None:
    metric_rows = [row for row in rows if row["metric"] == metric and row["status"] == "ok"]
    ordered = sorted(metric_rows, key=lambda row: float(row["p_raw"]))
    m = len(ordered)
    running_max = 0.0
    for idx, row in enumerate(ordered):
        adjusted = min(1.0, (m - idx) * float(row["p_raw"]))
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
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def median(values: list[float]) -> float:
    if not values:
        return float("nan")
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[midpoint]
    return (ordered[midpoint - 1] + ordered[midpoint]) / 2.0


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_latex_table(path: Path, rows: list[dict[str, object]], unit: str) -> None:
    ok_rows = [row for row in rows if row["status"] == "ok"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\\begin{table}[htbp]\n")
        handle.write("\\centering\n")
        if unit == "instance-mean":
            caption = (
                "Wilcoxon signed-rank test comparing CS-CGA with the compared methods "
                "on F1 score and SHD based on matched dataset-sample settings."
            )
        else:
            caption = (
                "Wilcoxon signed-rank test comparing CS-CGA with the compared methods "
                "on F1 score and SHD based on matched repeated runs."
            )
        handle.write(f"\\caption{{{caption}}}\n")
        handle.write("\\label{tab:wilcoxon_test}\n")
        handle.write("\\begin{tabular}{llrrrrr}\n")
        handle.write("\\toprule\n")
        handle.write("Compared method & Metric & Pairs & Mean diff. & $p$ & Holm-$p$ & Sig. \\\\\n")
        handle.write("\\midrule\n")
        for row in ok_rows:
            handle.write(
                f"{row['algorithm']} & {row['metric']} & {row['n_pairs']} & "
                f"{float(row['mean_improvement']):.4f} & {format_p(float(row['p_raw']))} & "
                f"{format_p(float(row['p_holm']))} & {row['significance']} \\\\\n"
            )
        handle.write("\\bottomrule\n")
        handle.write("\\end{tabular}\n")
        handle.write("\\\\\n")
        handle.write("\\footnotesize{Positive mean differences indicate better CS-CGA performance, i.e., higher F1 score or lower SHD.}\n")
        handle.write("\\end{table}\n")


def main() -> None:
    args = parse_args()
    result_dir = args.result_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    reference_index = build_file_index(result_dir / REFERENCE_ALGO)
    algorithm_indices = {
        algo_name: build_file_index(result_dir / dir_name)
        for algo_name, dir_name in COMPARISON_ALGOS.items()
    }

    summary_rows: list[dict[str, object]] = []
    detail_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []

    for algo_name, algo_index in algorithm_indices.items():
        metric_improvements = {"F1": [], "SHD": []}
        metric_reference_values = {"F1": [], "SHD": []}
        metric_comparison_values = {"F1": [], "SHD": []}
        matched_instances = 0

        for instance in sorted(reference_index.keys(), key=lambda item: (DATASETS.index(item[0]), item[1])):
            ref_path = reference_index[instance]
            cmp_path = algo_index.get(instance)
            dataset, sample_size = instance
            if cmp_path is None:
                missing_rows.append(
                    {
                        "algorithm": algo_name,
                        "dataset": dataset,
                        "sample_size": sample_size,
                        "reason": "missing_comparison_file",
                    }
                )
                continue

            ref_values = read_f1_shd(ref_path)
            cmp_values = read_f1_shd(cmp_path)
            n_pairs = min(len(ref_values["F1"]), len(cmp_values["F1"]), len(ref_values["SHD"]), len(cmp_values["SHD"]))
            if n_pairs < 2:
                missing_rows.append(
                    {
                        "algorithm": algo_name,
                        "dataset": dataset,
                        "sample_size": sample_size,
                        "reason": f"insufficient_pairs_{n_pairs}",
                    }
                )
                continue

            matched_instances += 1
            ref_f1 = ref_values["F1"][:n_pairs]
            cmp_f1 = cmp_values["F1"][:n_pairs]
            ref_shd = ref_values["SHD"][:n_pairs]
            cmp_shd = cmp_values["SHD"][:n_pairs]

            f1_improvements = [a - b for a, b in zip(ref_f1, cmp_f1)]
            shd_improvements = [b - a for a, b in zip(ref_shd, cmp_shd)]
            if args.unit == "run":
                metric_improvements["F1"].extend(f1_improvements)
                metric_improvements["SHD"].extend(shd_improvements)
                metric_reference_values["F1"].extend(ref_f1)
                metric_reference_values["SHD"].extend(ref_shd)
                metric_comparison_values["F1"].extend(cmp_f1)
                metric_comparison_values["SHD"].extend(cmp_shd)
            else:
                metric_improvements["F1"].append(mean(f1_improvements))
                metric_improvements["SHD"].append(mean(shd_improvements))
                metric_reference_values["F1"].append(mean(ref_f1))
                metric_reference_values["SHD"].append(mean(ref_shd))
                metric_comparison_values["F1"].append(mean(cmp_f1))
                metric_comparison_values["SHD"].append(mean(cmp_shd))

            detail_rows.append(
                {
                    "algorithm": algo_name,
                    "dataset": dataset,
                    "sample_size": sample_size,
                    "n_pairs": n_pairs,
                    "cscga_f1_mean": mean(ref_f1),
                    "comparison_f1_mean": mean(cmp_f1),
                    "f1_mean_improvement": mean(f1_improvements),
                    "cscga_shd_mean": mean(ref_shd),
                    "comparison_shd_mean": mean(cmp_shd),
                    "shd_mean_improvement": mean(shd_improvements),
                    "reference_file": str(ref_path),
                    "comparison_file": str(cmp_path),
                }
            )

        for metric in ["F1", "SHD"]:
            improvements = metric_improvements[metric]
            if len(improvements) < 2:
                summary_rows.append(
                    {
                        "algorithm": algo_name,
                        "metric": metric,
                        "status": "insufficient_data",
                        "analysis_unit": args.unit,
                        "matched_instances": matched_instances,
                        "n_pairs": len(improvements),
                    }
                )
                continue

            test = wilcoxon_signed_rank(improvements)
            wins = sum(1 for value in improvements if value > 1e-12)
            ties = sum(1 for value in improvements if abs(value) <= 1e-12)
            losses = sum(1 for value in improvements if value < -1e-12)
            summary_rows.append(
                {
                    "algorithm": algo_name,
                    "metric": metric,
                    "status": "ok",
                    "analysis_unit": args.unit,
                    "matched_instances": matched_instances,
                    "n_pairs": len(improvements),
                    "n_nonzero": test["n_nonzero"],
                    "cscga_mean": mean(metric_reference_values[metric]),
                    "comparison_mean": mean(metric_comparison_values[metric]),
                    "mean_improvement": mean(improvements),
                    "median_improvement": median(improvements),
                    "wins": wins,
                    "ties": ties,
                    "losses": losses,
                    "w_plus": test["w_plus"],
                    "w_minus": test["w_minus"],
                    "p_raw": test["p_raw"],
                    "p_method": test["p_method"],
                }
            )

    for metric in ["F1", "SHD"]:
        holm_adjust(summary_rows, metric)

    for row in summary_rows:
        if row.get("status") != "ok":
            continue
        p_holm = float(row["p_holm"])
        row["significance"] = significance_label(p_holm, args.alpha)
        row["interpretation"] = (
            "CS-CGA better" if float(row["mean_improvement"]) > 0 else "CS-CGA worse"
        )

    summary_fields = [
        "algorithm",
        "metric",
        "status",
        "analysis_unit",
        "matched_instances",
        "n_pairs",
        "n_nonzero",
        "cscga_mean",
        "comparison_mean",
        "mean_improvement",
        "median_improvement",
        "wins",
        "ties",
        "losses",
        "w_plus",
        "w_minus",
        "p_raw",
        "p_holm",
        "significance",
        "interpretation",
        "p_method",
    ]
    detail_fields = [
        "algorithm",
        "dataset",
        "sample_size",
        "n_pairs",
        "cscga_f1_mean",
        "comparison_f1_mean",
        "f1_mean_improvement",
        "cscga_shd_mean",
        "comparison_shd_mean",
        "shd_mean_improvement",
        "reference_file",
        "comparison_file",
    ]
    missing_fields = ["algorithm", "dataset", "sample_size", "reason"]

    write_csv(output_dir / "wilcoxon_paper_summary.csv", summary_rows, summary_fields)
    write_csv(output_dir / "wilcoxon_paper_instance_details.csv", detail_rows, detail_fields)
    write_csv(output_dir / "wilcoxon_paper_missing.csv", missing_rows, missing_fields)
    write_latex_table(output_dir / "wilcoxon_paper_table.tex", summary_rows, args.unit)

    print(f"Wrote: {output_dir / 'wilcoxon_paper_summary.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_paper_instance_details.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_paper_missing.csv'}")
    print(f"Wrote: {output_dir / 'wilcoxon_paper_table.tex'}")


if __name__ == "__main__":
    main()
