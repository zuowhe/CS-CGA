"""Main runner for pgmpy-based GTT, BS, and GES baselines.

The output CSV format follows the existing MATLAB baseline files:
trial, F1, sensitivity, specificity, precision, SHD, score, elapsed time,
steps, TP, FN, FP, TN, score_chushi.

Edit DEFAULT_INPUTS below for PyCharm runs, or pass dataset .mat files on the
command line.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

warnings.filterwarnings("ignore", message="`pgmpy.estimators.StructureScore` is deprecated.*", category=FutureWarning)
warnings.filterwarnings("ignore", message="GES is deprecated.*", category=FutureWarning)

import numpy as np
import pandas as pd
from pgmpy.estimators import BIC, GES

from open_source_bic_baselines import arcs_to_adjacency, load_frame, write_outputs
from paper_style_bic_baselines import (
    Timeout,
    bs_search,
    eval_against,
    gtt_search,
    normalize_network_name,
    parse_network,
)
from pgmpy_matlab_style_baselines import PgmpyBICScorer


ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = ROOT / "BS_GTT_standalone_package_20260618" / "datasets"

NETWORK_ORDER = [
    "Asia",
    "INS",
    "Water",
    "Alarm",
    "Hailfinder",
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


# PyCharm-friendly default requested by the manuscript experiments.
DEFAULT_INPUTS = ordered_default_inputs()

DEFAULT_ALGORITHMS = ["gtt", "bs"]
# DEFAULT_ALGORITHMS = ["gtt", "bs", "ges"]
DEFAULT_TRIALS = 10
DEFAULT_MAX_INDEGREE = 7
DEFAULT_PRUNE_MAX_PARENTS = 7
LARGE_TIMEOUT_NETWORKS = {"alarm", "hailfinder", "hepar", "win95pts", "and", "andes"}


@dataclass
class BaselineResult:
    algorithm: str
    dag: np.ndarray
    score: float
    elapsed_seconds: float
    steps: float
    timeout: bool
    extra: dict


def run_gtt(frame: pd.DataFrame, node_sizes: np.ndarray, args: argparse.Namespace) -> BaselineResult:
    scorer = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache)
    timeout = Timeout(args.timeout_seconds)
    result = gtt_search(scorer, timeout, args.max_indegree, args.min_improvement)
    return BaselineResult(
        algorithm="gtt",
        dag=result.dag,
        score=result.score,
        elapsed_seconds=result.elapsed_seconds,
        steps=float(result.steps),
        timeout=result.timeout,
        extra=result.extra,
    )


def run_bs(frame: pd.DataFrame, node_sizes: np.ndarray, args: argparse.Namespace) -> BaselineResult:
    scorer = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache)
    timeout = Timeout(args.timeout_seconds)
    result = bs_search(
        scorer,
        timeout,
        args.max_indegree,
        args.restarts,
        args.perturbations,
        args.seed,
        args.min_improvement,
    )
    return BaselineResult(
        algorithm="bs",
        dag=result.dag,
        score=result.score,
        elapsed_seconds=result.elapsed_seconds,
        steps=float(result.steps),
        timeout=result.timeout,
        extra=result.extra,
    )


def run_ges(frame: pd.DataFrame, node_sizes: np.ndarray, args: argparse.Namespace) -> BaselineResult:
    names = [str(column) for column in frame.columns]
    state_names = {
        names[i]: list(range(1, int(node_sizes[i]) + 1))
        for i in range(len(names))
    }
    scorer = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache)
    bic = BIC(frame, state_names=state_names)
    estimator = GES(frame, use_cache=args.use_cache, state_names=state_names)

    started = time.perf_counter()
    pdag = estimator.estimate(scoring_method=bic, min_improvement=args.ges_min_improvement, debug=args.ges_debug)
    dag_model = pdag.to_dag() if hasattr(pdag, "to_dag") else pdag
    elapsed = time.perf_counter() - started

    dag = arcs_to_adjacency(dag_model.edges(), names)
    return BaselineResult(
        algorithm="ges",
        dag=dag,
        score=scorer.dag_score(dag),
        elapsed_seconds=elapsed,
        steps=float("nan"),
        timeout=False,
        extra={
            "min_improvement": args.ges_min_improvement,
            "pdag_edges": int(len(list(pdag.edges()))),
            "dag_extension_edges": int(dag.sum()),
        },
    )


RUNNERS = {
    "gtt": run_gtt,
    "bs": run_bs,
    "ges": run_ges,
}

DISPLAY_NAME = {
    "gtt": "GTT_BIC",
    "bs": "BayesianSearch",
    "ges": "GES_BIC",
}


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
        f"{row[13]:11.3f}",
    ]


def print_matlab_header(dataset: str, algorithm: str) -> None:
    print(f"Running... {DISPLAY_NAME[algorithm]}_{dataset}", flush=True)
    print("Iter   F1 Score   Sensitivity   Specificity   Precision   HD    Bayes Score      Exe Time  #Steps TP  FN  FP  TN   score_chushi", flush=True)


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
            f"{int(row[12]):4d}   {row[13]:11.3f}",
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
    avg_row = ["Avg", *avg, float("nan")]
    std_row = ["Std", *std, float("nan")]
    return avg_row, std_row


def output_path(args: argparse.Namespace, dataset: str, algorithm: str) -> Path:
    folder = args.output_dir / DISPLAY_NAME[algorithm]
    folder.mkdir(parents=True, exist_ok=True)
    return folder / f"{dataset}_100_200_{DISPLAY_NAME[algorithm]}.csv"


def dataset_network_key(dataset: str) -> str:
    return re.sub(r"\d+$", "", dataset).lower()


def timeout_for_dataset(dataset: str, args: argparse.Namespace) -> float:
    network = dataset_network_key(dataset)
    if network in LARGE_TIMEOUT_NETWORKS:
        return args.large_timeout_seconds
    return args.timeout_seconds


def stable_prune_seed(dataset: str, algorithm: str, trial: int, base_seed: int) -> int:
    text = f"{dataset}:{algorithm}:{trial}"
    return int(base_seed + sum((idx + 1) * ord(ch) for idx, ch in enumerate(text)))


def random_prune_excess_parents(dag: np.ndarray, max_parents: int | None, seed: int) -> tuple[np.ndarray, int]:
    if max_parents is None or max_parents <= 0:
        return dag, 0
    pruned = dag.copy()
    rng = np.random.default_rng(seed)
    removed = 0
    for child in range(pruned.shape[1]):
        parents = np.flatnonzero(pruned[:, child])
        excess = len(parents) - max_parents
        if excess > 0:
            remove = rng.choice(parents, size=excess, replace=False)
            pruned[remove, child] = 0
            removed += int(excess)
    return pruned, removed


def run_dataset_algorithm(input_path: Path, algorithm: str, args: argparse.Namespace) -> None:
    rows = []
    dataset = input_path.stem
    dataset_timeout = timeout_for_dataset(dataset, args)
    print_matlab_header(dataset, algorithm)
    print(
        f"[config] timeout={dataset_timeout:.0f}s cache={'on' if args.use_cache else 'off'} "
        f"min_improvement={args.min_improvement} search_max_indegree={args.max_indegree} "
        f"random_prune_max_parents={args.prune_max_parents}",
        flush=True,
    )
    csv_path = output_path(args, dataset, algorithm)
    csv_path.unlink(missing_ok=True)

    for trial in range(1, args.trials + 1):
        print(f"[start] {DISPLAY_NAME[algorithm]} {dataset} trial {trial}/{args.trials}", flush=True)
        frame, _meta = load_frame(input_path, args.mat_var, trial, False)
        network_name = args.network or normalize_network_name(args.mat_var or input_path.stem)
        data_for_sizes = np.asarray(frame.astype(int).to_numpy().T, dtype=np.int64)
        truth, node_sizes = parse_network(network_name, args.root, data_for_sizes)

        run_args = argparse.Namespace(**vars(args))
        run_args.timeout_seconds = dataset_timeout
        result = RUNNERS[algorithm](frame, node_sizes, run_args)
        if result.timeout:
            timeout_row = [trial, "Timeout", result.elapsed_seconds, dataset_timeout]
            with csv_path.open("a", newline="") as handle:
                csv.writer(handle).writerow(timeout_row)
            print(f"{trial:4d}  Timeout after {result.elapsed_seconds:.3f} seconds", flush=True)
            break

        pruned_dag, removed_edges = random_prune_excess_parents(
            result.dag,
            args.prune_max_parents,
            stable_prune_seed(dataset, algorithm, trial, args.seed),
        )
        if removed_edges:
            score_started = time.perf_counter()
            result.dag = pruned_dag
            result.score = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache).dag_score(pruned_dag)
            result.elapsed_seconds += time.perf_counter() - score_started
            result.extra["random_parent_prune_removed"] = removed_edges
            print(f"[prune] removed {removed_edges} excess parent edge(s), MP={args.prune_max_parents}", flush=True)

        metrics = eval_against(result.dag, truth)
        row = matlab_style_line(trial, metrics, result.score, result.elapsed_seconds, result.steps)
        rows.append(row)
        print_matlab_row(row)

        names = [str(column) for column in frame.columns]
        edges = [(names[parent], names[child]) for parent, child in zip(*np.nonzero(result.dag))]
        write_outputs(args.graph_dir, f"{dataset}_trial{trial}", f"pgmpy_{algorithm}", names, edges, result.dag)

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
    parser.add_argument("--algorithms", nargs="+", choices=["gtt", "bs", "ges"], default=DEFAULT_ALGORITHMS)
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS)
    parser.add_argument("--mat-var", default=None)
    parser.add_argument("--network", default=None)
    parser.add_argument("--timeout-seconds", type=float, default=600, help="Timeout for smaller networks in seconds.")
    parser.add_argument("--large-timeout-seconds", type=float, default=7200, help="Timeout for Alarm and larger networks in seconds.")
    parser.add_argument("--max-indegree", type=int, default=DEFAULT_MAX_INDEGREE, help="Parent limit enforced during deletion/pruning.")
    parser.add_argument("--prune-max-parents", type=int, default=DEFAULT_PRUNE_MAX_PARENTS, help="Randomly delete excess parents after search. Use 0 to disable.")
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--perturbations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-improvement", type=float, default=1e-6)
    parser.add_argument("--use-cache", action="store_true", help="Use pgmpy ScoreCache. Default is off to better mimic slow MATLAB-style scans.")
    parser.add_argument("--ges-min-improvement", type=float, default=1e-6)
    parser.add_argument("--ges-debug", action="store_true")
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
        for algorithm in args.algorithms:
            run_dataset_algorithm(input_path, algorithm, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
