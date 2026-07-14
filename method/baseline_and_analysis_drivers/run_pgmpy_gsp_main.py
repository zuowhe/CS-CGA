"""Main runner for a pgmpy-scored Greedy Sparsest Permutation baseline."""

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
warnings.filterwarnings("ignore", message="`g_sq` is deprecated.*", category=FutureWarning)
warnings.filterwarnings("ignore", message="`power_divergence` is deprecated.*", category=FutureWarning)

import numpy as np
from scipy.stats import chi2
from pgmpy.estimators import CITests

from open_source_bic_baselines import load_frame, write_outputs
from paper_style_bic_baselines import normalize_network_name, parse_network
from pgmpy_matlab_style_baselines import PgmpyBICScorer


ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = ROOT / "GSP_standalone_package_20260618" / "datasets"
if not DATASET_DIR.exists():
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
LARGE_TIMEOUT_NETWORKS = {"alarm", "hailfinder", "hepar", "win95pts", "and", "andes"}

DISPLAY_NAME = "GSP"
DEFAULT_TRIALS = 10
DEFAULT_ALPHA = 0.01
DEFAULT_NUM_RUNS = 10
DEFAULT_SEARCH_DEPTH = 4


@dataclass
class GSPResult:
    dag: np.ndarray
    score: float
    elapsed_seconds: float
    steps: int
    timeout: bool
    extra: dict


class RuntimeLimit:
    def __init__(self, seconds: float):
        self.seconds = float(seconds)
        self.started = time.perf_counter()

    def expired(self) -> bool:
        return math.isfinite(self.seconds) and (time.perf_counter() - self.started) >= self.seconds

    def elapsed(self) -> float:
        return time.perf_counter() - self.started


class DiscreteLRTCI:
    def __init__(self, data: np.ndarray, node_sizes: np.ndarray, alpha: float, use_cache: bool = True):
        self.data = np.asarray(data, dtype=np.int64)
        self.node_sizes = np.asarray(node_sizes, dtype=np.int64)
        self.alpha = float(alpha)
        self.use_cache = use_cache
        self.cache: dict[tuple[int, int, tuple[int, ...]], bool] = {}
        self.requests = 0
        self.computed = 0

    def is_independent(self, x: int, y: int, cond_set: list[int] | tuple[int, ...] | np.ndarray) -> bool:
        x = int(x)
        y = int(y)
        z = tuple(sorted(int(node) for node in cond_set if int(node) not in (x, y)))
        key = (min(x, y), max(x, y), z)
        self.requests += 1
        if self.use_cache and key in self.cache:
            return self.cache[key]
        result = self._compute_lrt(x, y, z)
        if self.use_cache:
            self.cache[key] = result
        self.computed += 1
        return result

    def _compute_lrt(self, x: int, y: int, z: tuple[int, ...]) -> bool:
        sample_count = self.data.shape[0]
        x_states = int(self.node_sizes[x])
        y_states = int(self.node_sizes[y])
        base_df = (x_states - 1) * (y_states - 1)
        z_configs = 1
        for node in z:
            z_configs *= int(self.node_sizes[node])
            if base_df * z_configs > sample_count:
                return False
        df = base_df * z_configs
        if df <= 0 or sample_count < df:
            return False

        x_values = self.data[:, x] - 1
        y_values = self.data[:, y] - 1
        if not z:
            table = np.zeros((x_states, y_states), dtype=float)
            valid = (
                (x_values >= 0)
                & (x_values < x_states)
                & (y_values >= 0)
                & (y_values < y_states)
            )
            np.add.at(table, (x_values[valid], y_values[valid]), 1.0)
            statistic = lrt_statistic(table)
        else:
            z_values = self.data[:, z] - 1
            valid = (
                (x_values >= 0)
                & (x_values < x_states)
                & (y_values >= 0)
                & (y_values < y_states)
                & np.all(z_values >= 0, axis=1)
            )
            for pos, node in enumerate(z):
                valid &= z_values[:, pos] < int(self.node_sizes[node])
            if not np.any(valid):
                return False
            _, group_id = np.unique(z_values[valid, :], axis=0, return_inverse=True)
            valid_x = x_values[valid]
            valid_y = y_values[valid]
            statistic = 0.0
            for group in range(int(group_id.max()) + 1):
                rows = group_id == group
                table = np.zeros((x_states, y_states), dtype=float)
                np.add.at(table, (valid_x[rows], valid_y[rows]), 1.0)
                statistic += lrt_statistic(table)

        if not math.isfinite(statistic):
            return False
        p_value = float(chi2.sf(statistic, df))
        return p_value >= self.alpha


class PgmpyGSqCI:
    def __init__(self, frame, alpha: float, use_cache: bool = True):
        self.frame = frame
        self.names = [str(column) for column in frame.columns]
        self.alpha = float(alpha)
        self.use_cache = use_cache
        self.cache: dict[tuple[int, int, tuple[int, ...]], bool] = {}
        self.requests = 0
        self.computed = 0

    def is_independent(self, x: int, y: int, cond_set: list[int] | tuple[int, ...] | np.ndarray) -> bool:
        x = int(x)
        y = int(y)
        z = tuple(sorted(int(node) for node in cond_set if int(node) not in (x, y)))
        key = (min(x, y), max(x, y), z)
        self.requests += 1
        if self.use_cache and key in self.cache:
            return self.cache[key]
        try:
            result = bool(
                CITests.g_sq(
                    X=self.names[x],
                    Y=self.names[y],
                    Z=[self.names[node] for node in z],
                    data=self.frame,
                    boolean=True,
                    significance_level=self.alpha,
                )
            )
        except Exception:
            result = False
        if self.use_cache:
            self.cache[key] = result
        self.computed += 1
        return result


def lrt_statistic(table: np.ndarray) -> float:
    total = float(table.sum())
    if total <= 0:
        return 0.0
    expected = np.outer(table.sum(axis=1), table.sum(axis=0)) / total
    mask = (table > 0) & (expected > 0)
    if not np.any(mask):
        return 0.0
    return float(2.0 * np.sum(table[mask] * np.log(table[mask] / expected[mask])))


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


def permutation_to_dag(perm: np.ndarray, ci: DiscreteLRTCI, limit: RuntimeLimit) -> tuple[np.ndarray, bool]:
    node_count = len(perm)
    dag = np.zeros((node_count, node_count), dtype=np.int8)
    for second_pos in range(node_count):
        if limit.expired():
            return dag, True
        child = int(perm[second_pos])
        predecessors = [int(node) for node in perm[:second_pos]]
        for first_pos in range(second_pos):
            parent = int(perm[first_pos])
            cond_set = [node for node in predecessors if node != parent]
            if not ci.is_independent(parent, child, cond_set):
                dag[parent, child] = 1
            if limit.expired():
                return dag, True
    return dag, False


def covered_arcs(dag: np.ndarray) -> list[tuple[int, int]]:
    arcs = []
    parents, children = np.nonzero(dag)
    for parent, child in zip(parents, children):
        parent_parents = set(int(node) for node in np.flatnonzero(dag[:, parent]))
        child_parents = set(int(node) for node in np.flatnonzero(dag[:, child]))
        child_parents.discard(int(parent))
        if parent_parents == child_parents:
            arcs.append((int(parent), int(child)))
    return arcs


def update_minimal_imap(dag: np.ndarray, parent: int, child: int, ci: DiscreteLRTCI, limit: RuntimeLimit) -> tuple[list[tuple[int, int]], bool]:
    parent_set = [int(node) for node in np.flatnonzero(dag[:, parent])]
    remove_arcs: set[tuple[int, int]] = set()
    for p_node in parent_set:
        if limit.expired():
            return sorted(remove_arcs), True
        cond_set = [node for node in parent_set if node != p_node] + [child]
        if ci.is_independent(parent, p_node, cond_set):
            remove_arcs.add((p_node, parent))

        if dag[p_node, child] != 0:
            cond_set = [node for node in parent_set if node != p_node]
            if ci.is_independent(child, p_node, cond_set):
                remove_arcs.add((p_node, child))
    return sorted(remove_arcs), False


def candidate_moves(dag: np.ndarray, ci: DiscreteLRTCI, limit: RuntimeLimit) -> tuple[list[dict], bool]:
    moves = []
    for parent, child in covered_arcs(dag):
        if limit.expired():
            return moves, True
        remove_arcs, timeout = update_minimal_imap(dag, parent, child, ci, limit)
        if timeout:
            return moves, True
        moves.append(
            {
                "parent": parent,
                "child": child,
                "remove_arcs": remove_arcs,
                "remove_count": len(remove_arcs),
            }
        )
    moves.sort(key=lambda move: move["remove_count"], reverse=True)
    return moves, False


def apply_gsp_move(dag: np.ndarray, move: dict) -> np.ndarray:
    candidate = dag.copy()
    parent = move["parent"]
    child = move["child"]
    candidate[parent, child] = 0
    candidate[child, parent] = 1
    for edge_parent, edge_child in move["remove_arcs"]:
        candidate[edge_parent, edge_child] = 0
    return candidate


def dag_key(dag: np.ndarray) -> bytes:
    return dag.astype(np.int8, copy=False).tobytes()


def depth_limited_sparser_search(root_dag: np.ndarray, ci: DiscreteLRTCI, limit: RuntimeLimit, search_depth: int) -> tuple[np.ndarray, bool, bool, int]:
    if search_depth <= 0:
        return root_dag, False, False, 0
    visited = {dag_key(root_dag): 0}
    return dfs_sparser(root_dag, 0, int(root_dag.sum()), ci, limit, search_depth, visited)


def dfs_sparser(
    current_dag: np.ndarray,
    depth: int,
    root_edges: int,
    ci: DiscreteLRTCI,
    limit: RuntimeLimit,
    search_depth: int,
    visited: dict[bytes, int],
) -> tuple[np.ndarray, bool, bool, int]:
    if depth >= search_depth:
        return current_dag, False, False, 0
    moves, timeout = candidate_moves(current_dag, ci, limit)
    if timeout:
        return current_dag, False, True, 0
    current_edges = int(current_dag.sum())
    for move in moves:
        if limit.expired():
            return current_dag, False, True, 0
        candidate = apply_gsp_move(current_dag, move)
        candidate_depth = depth + 1
        key = dag_key(candidate)
        if key in visited and visited[key] <= candidate_depth:
            continue
        visited[key] = candidate_depth
        candidate_edges = int(candidate.sum())
        if candidate_edges < root_edges:
            return candidate, True, False, 1
        if candidate_edges <= current_edges:
            found_dag, found, timeout, descendant_path = dfs_sparser(
                candidate, candidate_depth, root_edges, ci, limit, search_depth, visited
            )
            if timeout:
                return current_dag, False, True, 0
            if found:
                return found_dag, True, False, descendant_path + 1
    return current_dag, False, False, 0


def gsp_search(frame, data: np.ndarray, node_sizes: np.ndarray, scorer: PgmpyBICScorer, args: argparse.Namespace, rng: np.random.Generator) -> GSPResult:
    limit = RuntimeLimit(args.timeout_seconds)
    if args.ci_backend == "pgmpy_gsq":
        ci = PgmpyGSqCI(frame, args.alpha, use_cache=args.use_ci_cache)
    else:
        ci = DiscreteLRTCI(data, node_sizes, args.alpha, use_cache=args.use_ci_cache)
    node_count = len(node_sizes)
    terminal: list[tuple[np.ndarray, int]] = []
    steps = 0

    for _run in range(args.num_runs):
        if limit.expired():
            break
        perm = rng.permutation(node_count)
        current_dag, timeout = permutation_to_dag(perm, ci, limit)
        if timeout:
            break

        while not limit.expired():
            next_dag, found, timeout, path_length = depth_limited_sparser_search(
                current_dag, ci, limit, args.search_depth
            )
            if timeout:
                break
            if found:
                current_dag = next_dag
                steps += path_length
            else:
                terminal.append((current_dag, int(current_dag.sum())))
                break

    if not terminal:
        return GSPResult(
            dag=np.zeros((node_count, node_count), dtype=np.int8),
            score=float("-inf"),
            elapsed_seconds=limit.elapsed(),
            steps=steps,
            timeout=True,
            extra={"ci_requests": ci.requests, "ci_computed": ci.computed},
        )

    best_dag = None
    best_score = float("-inf")
    best_edges = math.inf
    for dag, edge_count in terminal:
        score = scorer.dag_score(dag)
        if edge_count < best_edges or (edge_count == best_edges and score > best_score):
            best_dag = dag
            best_score = score
            best_edges = edge_count

    return GSPResult(
        dag=best_dag,
        score=best_score,
        elapsed_seconds=limit.elapsed(),
        steps=steps,
        timeout=limit.expired(),
        extra={
            "ci_requests": ci.requests,
            "ci_computed": ci.computed,
            "edge_count": int(best_edges),
            "runs": len(terminal),
            "search_depth": args.search_depth,
        },
    )


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
            f"{row[8]:7.3f}",
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
        f"{int(row[8]):4d}",
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
        print(
            f"{int(row[0]):4d}  {row[1]:9.5f}    {row[2]:9.5f}    {row[3]:9.5f}  "
            f"{row[4]:9.5f}  {row[5]:3.1f}   {row[6]:11.3f}  {row[7]:11.3f}      "
            f"{int(row[8]):4d}       {int(row[9]):4d}{int(row[10]):4d}{int(row[11]):4d}"
            f"{int(row[12]):4d}{int(row[13]):4d}   {row[14]:11.3f}",
            flush=True,
        )


def summarize_rows(rows: list[list[float]]) -> tuple[list, list]:
    numeric = np.asarray([[float(value) for value in row[1:9]] for row in rows], dtype=float)
    avg = np.nanmean(numeric, axis=0)
    std = np.nanstd(numeric, axis=0, ddof=0)
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


def stable_seed(dataset: str, trial: int, base_seed: int) -> int:
    text = f"{dataset}:{trial}:gsp"
    return int(base_seed + sum((idx + 1) * ord(ch) for idx, ch in enumerate(text)))


def run_dataset(input_path: Path, args: argparse.Namespace) -> None:
    rows = []
    dataset = input_path.stem
    dataset_timeout = timeout_for_dataset(dataset, args)
    print_matlab_header(dataset)
    print(
        f"[config] timeout={dataset_timeout:.0f}s alpha={args.alpha:g} "
        f"num_runs={args.num_runs} search_depth={args.search_depth} "
        f"ci_backend={args.ci_backend} ci_cache={'on' if args.use_ci_cache else 'off'} "
        f"score_cache={'on' if args.use_cache else 'off'}",
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
        data = np.asarray(frame.astype(int).to_numpy(), dtype=np.int64)

        run_args = argparse.Namespace(**vars(args))
        run_args.timeout_seconds = dataset_timeout
        scorer = PgmpyBICScorer(frame, node_sizes, use_cache=args.use_cache)
        rng = np.random.default_rng(stable_seed(dataset, trial, args.seed))
        result = gsp_search(frame, data, node_sizes, scorer, run_args, rng)

        if result.timeout and not math.isfinite(result.score):
            timeout_row = [trial, "Timeout", result.elapsed_seconds, dataset_timeout]
            with csv_path.open("a", newline="") as handle:
                csv.writer(handle).writerow(timeout_row)
            print(f"{trial:4d}  Timeout after {result.elapsed_seconds:.3f} seconds", flush=True)
            break

        metrics = eval_against_matlab(result.dag, truth)
        row = matlab_style_line(trial, metrics, result.score, result.elapsed_seconds, result.steps)
        rows.append(row)
        print_matlab_row(row)

        names = [str(column) for column in frame.columns]
        edges = [(names[parent], names[child]) for parent, child in zip(*np.nonzero(result.dag))]
        write_outputs(args.graph_dir, f"{dataset}_trial{trial}", "pgmpy_gsp", names, edges, result.dag)

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
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("--ci-backend", choices=["matlab_lrt", "pgmpy_gsq"], default="matlab_lrt")
    parser.add_argument("--num-runs", type=int, default=DEFAULT_NUM_RUNS)
    parser.add_argument("--search-depth", type=int, default=DEFAULT_SEARCH_DEPTH)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--use-ci-cache", action="store_true")
    parser.add_argument("--use-cache", action="store_true")
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "pgmpy_matlab_style")
    parser.add_argument("--graph-dir", type=Path, default=ROOT / "results" / "pgmpy_matlab_style" / "graphs")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.root = args.root.resolve()
    args.num_runs = max(1, int(math.ceil(args.num_runs)))
    args.search_depth = max(0, int(math.ceil(args.search_depth)))
    inputs = args.inputs or DEFAULT_INPUTS
    for input_path in inputs:
        input_path = input_path if input_path.is_absolute() else (ROOT / input_path)
        if not input_path.exists():
            raise FileNotFoundError(input_path)
        run_dataset(input_path, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
