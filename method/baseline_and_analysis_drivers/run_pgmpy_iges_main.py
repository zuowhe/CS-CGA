"""Main runner for the pgmpy-based IGES_RCI baseline."""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="`pgmpy.estimators.StructureScore` is deprecated.*", category=FutureWarning)
warnings.filterwarnings("ignore", message="GES is deprecated.*", category=FutureWarning)

import numpy as np
import pandas as pd

from open_source_bic_baselines import load_frame, write_outputs
from paper_style_bic_baselines import normalize_network_name, parse_network
from pgmpy_matlab_style_baselines import PgmpyBICScorer
from run_pgmpy_gtt_bs_ges_main import (
    DATASET_DIR,
    LARGE_TIMEOUT_NETWORKS,
    ROOT,
    random_prune_excess_parents,
    run_ges,
    stable_prune_seed,
)


NETWORK_ORDER = [
    # "Asia",
    # "INS",
    # "Water",
    # "Alarm",
    # "Hailfinder",
    "HEPAR",
    "Win95pts",
    "AND",
]
SAMPLE_SIZE_ORDER = [500, 1000, 3000]


def ordered_default_inputs() -> list[Path]:
    inputs = []
    for network in NETWORK_ORDER:
        for sample_size in SAMPLE_SIZE_ORDER:
            path = DATASET_DIR / f"{network}{sample_size}.mat"
            if path.exists():
                inputs.append(path)
    known = {path.name for path in inputs}
    inputs.extend(path for path in sorted(DATASET_DIR.glob("*.mat")) if path.name not in known)
    return inputs


DEFAULT_INPUTS = ordered_default_inputs()
DEFAULT_TRIALS = 10
DEFAULT_MAX_INDEGREE = 7
DEFAULT_PRUNE_MAX_PARENTS = 7
DISPLAY_NAME = "IGES_RCI"


def joint_entropy(frame_values: np.ndarray, nodes: list[int] | np.ndarray) -> float:
    nodes = sorted({int(node) for node in nodes})
    if not nodes:
        return 0.0
    observations = frame_values[:, nodes].astype(float, copy=False)
    finite = np.isfinite(observations).all(axis=1)
    observations = observations[finite]
    sample_count = observations.shape[0]
    if sample_count == 0:
        return 0.0
    _, counts = np.unique(observations, axis=0, return_counts=True)
    prob = counts[counts > 0] / sample_count
    return float(-np.sum(prob * np.log2(prob)))


def conditional_mutual_information(frame_values: np.ndarray, x: int, y: int, z_set: list[int] | np.ndarray) -> float:
    z_nodes = sorted({int(node) for node in z_set})
    if not z_nodes:
        cmi = joint_entropy(frame_values, [x]) + joint_entropy(frame_values, [y]) - joint_entropy(frame_values, [x, y])
    else:
        cmi = (
            joint_entropy(frame_values, [x, *z_nodes])
            + joint_entropy(frame_values, [y, *z_nodes])
            - joint_entropy(frame_values, z_nodes)
            - joint_entropy(frame_values, [x, y, *z_nodes])
        )
    if not math.isfinite(cmi) or cmi < 0:
        cmi = max(cmi, 0.0)
    return float(cmi)


def pairwise_mutual_information(frame_values: np.ndarray) -> np.ndarray:
    node_count = frame_values.shape[1]
    mi_pair = np.zeros((node_count, node_count), dtype=float)
    for i in range(node_count - 1):
        for j in range(i + 1, node_count):
            mi_pair[i, j] = conditional_mutual_information(frame_values, i, j, [])
            mi_pair[j, i] = mi_pair[i, j]
    return mi_pair


def iges_cmi_ib_filter(
    dag: np.ndarray,
    frame: pd.DataFrame,
    cmi_threshold: float | None,
    ib_loss_threshold: float | None,
    cmi_percentile: float | None,
    ib_loss_percentile: float | None,
    beta: float,
) -> tuple[np.ndarray, dict]:
    filtered = dag.copy()
    parents, children = np.nonzero(dag)
    edge_count = len(parents)
    info = {
        "removed_edges": 0,
        "cmi_removed": 0,
        "ib_removed": 0,
        "edge_count_before": int(edge_count),
        "edge_count_after": int(edge_count),
        "cmi_threshold": cmi_threshold,
        "ib_loss_threshold": ib_loss_threshold,
        "cmi_percentile": cmi_percentile,
        "ib_loss_percentile": ib_loss_percentile,
        "beta": beta,
    }
    if edge_count == 0:
        return filtered, info

    frame_values = frame.astype(float).to_numpy()
    mi_pair = pairwise_mutual_information(frame_values)
    cmi_scores = np.zeros(edge_count, dtype=float)
    ib_losses = np.zeros(edge_count, dtype=float)

    for edge_id, (parent, child) in enumerate(zip(parents, children)):
        cond_set = [int(node) for node in np.flatnonzero(dag[:, child]) if int(node) != int(parent)]
        cmi_scores[edge_id] = conditional_mutual_information(frame_values, int(parent), int(child), cond_set)
        redundancy = float(np.mean(mi_pair[int(parent), cond_set])) if cond_set else 0.0
        ib_losses[edge_id] = beta * redundancy / (cmi_scores[edge_id] + np.finfo(float).eps)

    if cmi_threshold is None:
        cmi_threshold = float(np.percentile(cmi_scores, cmi_percentile if cmi_percentile is not None else 10.0))
    if ib_loss_threshold is None:
        ib_loss_threshold = float(np.percentile(ib_losses, ib_loss_percentile if ib_loss_percentile is not None else 90.0))
    info["cmi_threshold"] = cmi_threshold
    info["ib_loss_threshold"] = ib_loss_threshold

    remove_by_cmi = cmi_scores < cmi_threshold
    remove_by_ib = ib_losses > ib_loss_threshold
    remove_flags = remove_by_cmi | remove_by_ib
    for parent, child, remove in zip(parents, children, remove_flags):
        if remove:
            filtered[int(parent), int(child)] = 0

    info["cmi_removed"] = int(np.count_nonzero(remove_by_cmi))
    info["ib_removed"] = int(np.count_nonzero(remove_by_ib))
    info["removed_edges"] = int(np.count_nonzero(remove_flags))
    info["edge_count_after"] = int(filtered.sum())
    return filtered, info


def eval_against_matlab(pred: np.ndarray, truth: np.ndarray | None) -> dict:
    if truth is None:
        return {}
    n = truth.shape[0]
    tp = tp2 = fn = fp = tn = 0
    for i in range(n):
        for j in range(i, n):
            t_ij, t_ji = truth[i, j], truth[j, i]
            p_ij, p_ji = pred[i, j], pred[j, i]
            if t_ij or t_ji:
                if (t_ij and p_ij) or (t_ji and p_ji):
                    tp += 1
                elif (t_ij and p_ji) or (t_ji and p_ij):
                    tp2 += 1
                else:
                    fn += 1
            else:
                if p_ij or p_ji:
                    fp += 1
                else:
                    tn += 1
    weighted_tp = tp + tp2 / 2
    sensitivity = weighted_tp / (weighted_tp + fn) if (weighted_tp + fn) else 1.0
    precision = weighted_tp / (weighted_tp + fp) if (weighted_tp + fp) else 1.0
    f1 = 2 * sensitivity * precision / (sensitivity + precision) if (sensitivity + precision) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    shd = fn + fp + tp2 / 2
    return {
        "f1": f1,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "shd": shd,
        "tp": tp,
        "tp2": tp2,
        "fn": fn,
        "fp": fp,
        "tn": tn,
    }


def run_iges(frame: pd.DataFrame, node_sizes: np.ndarray, args: argparse.Namespace):
    result = run_ges(frame, node_sizes, args)
    result.algorithm = "iges"
    if result.timeout:
        return result

    cmi_threshold = args.iges_cmi_threshold
    ib_loss_threshold = args.iges_ib_loss_threshold
    cmi_percentile = args.iges_cmi_percentile
    ib_loss_percentile = args.iges_ib_loss_percentile
    if args.iges_filter == "fixed":
        cmi_threshold = 1.0 if cmi_threshold is None else cmi_threshold
        ib_loss_threshold = 10000.0 if ib_loss_threshold is None else ib_loss_threshold
        cmi_percentile = None
        ib_loss_percentile = None

    if args.iges_filter != "none":
        started = time.perf_counter()
        filtered_dag, filter_info = iges_cmi_ib_filter(
            result.dag,
            frame,
            cmi_threshold,
            ib_loss_threshold,
            cmi_percentile,
            ib_loss_percentile,
            args.iges_beta,
        )
        scorer = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache)
        result.dag = filtered_dag
        result.score = scorer.dag_score(filtered_dag)
        result.elapsed_seconds += time.perf_counter() - started
        if not math.isnan(result.steps):
            result.steps += filter_info["removed_edges"]
        result.extra = {**result.extra, "iges_filter": {**filter_info, "mode": args.iges_filter}}
    else:
        result.extra = {**result.extra, "iges_filter": {"mode": "none", "removed_edges": 0}}
    return result


def matlab_style_line(trial, metrics, score, elapsed, steps, score_chushi=-1.0):
    return [
        trial,
        metrics["f1"],
        metrics["sensitivity"],
        metrics["specificity"],
        metrics["precision"],
        metrics["shd"],
        score,
        elapsed,
        steps,
        metrics["tp"],
        metrics["tp2"],
        metrics["fn"],
        metrics["fp"],
        metrics["tn"],
        score_chushi,
    ]


def format_csv_row(row) -> list[str]:
    label = row[0]
    if isinstance(label, str):
        return [
            label,
            f"{row[1]:9.5f}",
            f"{row[2]:9.5f}",
            f"{row[3]:9.5f}",
            f"{row[4]:9.5f}",
            f"{row[5]:5.3f}",
            f"{row[6]:11.3f}",
            f"{row[7]:11.3f}",
            f"{row[8]:7.3f}" if not math.isnan(float(row[8])) else "    NaN",
            "       NaN",
        ]
    return [
        f"{int(row[0]):4d}",
        f"{row[1]:9.5f}",
        f"{row[2]:9.5f}",
        f"{row[3]:9.5f}",
        f"{row[4]:9.5f}",
        f"{row[5]:3.1f}",
        f"{row[6]:11.3f}",
        f"{row[7]:11.3f}",
        f"{int(row[8]):4d}" if not math.isnan(float(row[8])) else " NaN",
        f"{int(row[9]):4d}",
        f"{int(row[10]):4d}",
        f"{int(row[11]):4d}",
        f"{int(row[12]):4d}",
        f"{int(row[13]):4d}",
        f"{row[14]:11.3f}",
    ]


def print_matlab_header(dataset: str) -> None:
    print(f"Running... {DISPLAY_NAME}_{dataset}", flush=True)
    print("Iter   F1 Score   Sensitivity   Specificity   Precision   HD    Bayes Score      Exe Time  #Steps TP  TP2  FN  FP  TN   score_chushi", flush=True)


def print_matlab_row(row) -> None:
    if isinstance(row[0], str):
        print(
            f"{row[0]:>4s}  {row[1]:9.5f}    {row[2]:9.5f}    {row[3]:9.5f}  "
            f"{row[4]:9.5f}  {row[5]:5.3f}   {row[6]:11.3f}  {row[7]:11.3f}      "
            f"{row[8]:7.3f}",
            flush=True,
        )
    else:
        steps = f"{int(row[8]):4d}" if not math.isnan(float(row[8])) else " NaN"
        print(
            f"{int(row[0]):4d}  {row[1]:9.5f}    {row[2]:9.5f}    {row[3]:9.5f}  "
            f"{row[4]:9.5f}  {row[5]:3.1f}   {row[6]:11.3f}  {row[7]:11.3f}      "
            f"{steps}       {int(row[9]):4d}{int(row[10]):4d}{int(row[11]):4d}"
            f"{int(row[12]):4d}{int(row[13]):4d}   {row[14]:11.3f}",
            flush=True,
        )


def summarize_rows(rows: list[list[float]]) -> tuple[list, list]:
    numeric = np.asarray([[float(value) for value in row[1:9]] for row in rows], dtype=float)
    avg = np.array([
        float("nan") if np.all(np.isnan(column)) else float(np.nanmean(column))
        for column in numeric.T
    ])
    std = np.array([
        float("nan") if np.all(np.isnan(column)) else float(np.nanstd(column, ddof=0))
        for column in numeric.T
    ])
    return ["Avg", *avg, float("nan")], ["Std", *std, float("nan")]


def output_path(args: argparse.Namespace, dataset: str) -> Path:
    folder = args.output_dir / DISPLAY_NAME
    folder.mkdir(parents=True, exist_ok=True)
    return folder / f"{dataset}_100_200_{DISPLAY_NAME}.csv"


def dataset_network_key(dataset: str) -> str:
    return re.sub(r"\d+$", "", dataset).lower()


def timeout_for_dataset(dataset: str, args: argparse.Namespace) -> float:
    network = dataset_network_key(dataset)
    if network in LARGE_TIMEOUT_NETWORKS:
        return args.large_timeout_seconds
    return args.timeout_seconds


def run_dataset(input_path: Path, args: argparse.Namespace) -> None:
    rows = []
    dataset = input_path.stem
    dataset_timeout = timeout_for_dataset(dataset, args)
    print_matlab_header(dataset)
    print(
        f"[config] timeout={dataset_timeout:.0f}s cache={'on' if args.use_cache else 'off'} "
        f"filter={args.iges_filter} cmi_threshold={args.iges_cmi_threshold} "
        f"ib_loss_threshold={args.iges_ib_loss_threshold} "
        f"cmi_percentile={args.iges_cmi_percentile:g} "
        f"ib_loss_percentile={args.iges_ib_loss_percentile:g} beta={args.iges_beta:g} "
        f"random_prune_max_parents={args.prune_max_parents}",
        flush=True,
    )
    csv_path = output_path(args, dataset)
    csv_path.unlink(missing_ok=True)

    for trial in range(1, args.trials + 1):
        print(f"[start] {DISPLAY_NAME} {dataset} trial {trial}/{args.trials}", flush=True)
        frame, _meta = load_frame(input_path, args.mat_var, trial, False)
        network_name = args.network or normalize_network_name(args.mat_var or input_path.stem)
        data_for_sizes = np.asarray(frame.astype(int).to_numpy().T, dtype=np.int64)
        truth, node_sizes = parse_network(network_name, args.root, data_for_sizes)

        run_args = argparse.Namespace(**vars(args))
        run_args.timeout_seconds = dataset_timeout
        result = run_iges(frame, node_sizes, run_args)
        if result.timeout:
            timeout_row = [trial, "Timeout", result.elapsed_seconds, dataset_timeout]
            with csv_path.open("a", newline="") as handle:
                csv.writer(handle).writerow(timeout_row)
            print(f"{trial:4d}  Timeout after {result.elapsed_seconds:.3f} seconds", flush=True)
            break

        pruned_dag, removed_edges = random_prune_excess_parents(
            result.dag,
            args.prune_max_parents,
            stable_prune_seed(dataset, "iges", trial, args.seed),
        )
        if removed_edges:
            score_started = time.perf_counter()
            result.dag = pruned_dag
            result.score = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache).dag_score(pruned_dag)
            result.elapsed_seconds += time.perf_counter() - score_started
            result.extra["random_parent_prune_removed"] = removed_edges
            print(f"[prune] removed {removed_edges} excess parent edge(s), MP={args.prune_max_parents}", flush=True)

        metrics = eval_against_matlab(result.dag, truth)
        row = matlab_style_line(trial, metrics, result.score, result.elapsed_seconds, result.steps)
        rows.append(row)
        print_matlab_row(row)

        names = [str(column) for column in frame.columns]
        edges = [(names[parent], names[child]) for parent, child in zip(*np.nonzero(result.dag))]
        write_outputs(args.graph_dir, f"{dataset}_trial{trial}", "pgmpy_iges", names, edges, result.dag)

        with csv_path.open("a", newline="") as handle:
            csv.writer(handle).writerow(format_csv_row(row))

    if rows:
        avg_row, std_row = summarize_rows(rows)
        print_matlab_row(avg_row)
        print_matlab_row(std_row)
        with csv_path.open("a", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(format_csv_row(avg_row))
            writer.writerow(format_csv_row(std_row))

    print(f"[saved] {csv_path}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="*", type=Path, help="Dataset .mat files. If omitted, DEFAULT_INPUTS in this script are used.")
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS)
    parser.add_argument("--mat-var", default=None)
    parser.add_argument("--network", default=None)
    parser.add_argument("--timeout-seconds", type=float, default=600)
    parser.add_argument("--large-timeout-seconds", type=float, default=7200)
    parser.add_argument("--max-indegree", type=int, default=DEFAULT_MAX_INDEGREE)
    parser.add_argument("--prune-max-parents", type=int, default=DEFAULT_PRUNE_MAX_PARENTS)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--use-cache", action="store_true")
    parser.add_argument("--ges-min-improvement", type=float, default=1e-6)
    parser.add_argument("--ges-debug", action="store_true")
    parser.add_argument("--iges-filter", choices=["none", "percentile", "fixed"], default="percentile")
    parser.add_argument("--iges-cmi-threshold", type=float, default=None)
    parser.add_argument("--iges-ib-loss-threshold", type=float, default=None)
    parser.add_argument("--iges-cmi-percentile", type=float, default=10.0)
    parser.add_argument("--iges-ib-loss-percentile", type=float, default=90.0)
    parser.add_argument("--iges-beta", type=float, default=1.0)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "pgmpy_matlab_style")
    parser.add_argument("--graph-dir", type=Path, default=ROOT / "results" / "pgmpy_matlab_style" / "graphs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.root = args.root.resolve()
    inputs = args.inputs or DEFAULT_INPUTS
    for input_path in inputs:
        input_path = input_path if input_path.is_absolute() else (ROOT / input_path)
        if not input_path.exists():
            raise FileNotFoundError(input_path)
        run_dataset(input_path, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
