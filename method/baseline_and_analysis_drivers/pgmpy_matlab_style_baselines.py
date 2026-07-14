"""Run MATLAB-style GTT/BS search using pgmpy local BIC scores.

This keeps pgmpy as the open-source scoring backend, but mirrors the search
procedures used by the MATLAB baselines:

- GTT: greedy thickening followed by greedy thinning.
- BS: maximum-branching initialization, local add/delete/reverse hill-climb,
  random perturbation, and restarts.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from pgmpy.estimators import BIC
from pgmpy.estimators.ScoreCache import ScoreCache

from open_source_bic_baselines import arcs_to_adjacency, load_frame, write_outputs
from paper_style_bic_baselines import (
    bs_search,
    eval_against,
    gtt_search,
    matlab_row,
    normalize_network_name,
    parse_network,
    Timeout,
)


class PgmpyBICScorer:
    def __init__(self, frame: pd.DataFrame, node_sizes: np.ndarray, use_cache: bool = True):
        self.frame = frame
        self.names = [str(column) for column in frame.columns]
        self.n_nodes = len(self.names)
        state_names = {
            self.names[i]: list(range(1, int(node_sizes[i]) + 1))
            for i in range(self.n_nodes)
        }
        base = BIC(frame, state_names=state_names)
        self.score_obj = ScoreCache(base, frame, max_size=200000, state_names=state_names) if use_cache else base

    def local_score(self, child: int, parents: list[int] | tuple[int, ...] | np.ndarray) -> float:
        parent_names = [self.names[int(parent)] for parent in parents]
        return float(self.score_obj.local_score(self.names[int(child)], parent_names))

    def dag_score(self, dag: np.ndarray) -> float:
        return float(sum(self.local_score(child, np.flatnonzero(dag[:, child])) for child in range(self.n_nodes)))


def strict_eval_against(dag: np.ndarray, truth: np.ndarray) -> dict:
    return eval_against(dag, truth)


def run_one(args: argparse.Namespace, input_path: Path, algorithm: str) -> dict:
    frame, meta = load_frame(input_path, args.mat_var, args.trial, args.csv_nodes_in_rows)
    compare_root = args.compare_root.resolve()
    network_name = args.network or normalize_network_name(args.mat_var or input_path.stem)
    data_for_sizes = np.asarray(frame.astype(int).to_numpy().T, dtype=np.int64)
    truth, node_sizes = parse_network(network_name, compare_root, data_for_sizes)
    scorer = PgmpyBICScorer(frame, node_sizes, use_cache=not args.no_cache)
    timeout = Timeout(args.timeout_seconds)

    if algorithm == "gtt":
        result = gtt_search(scorer, timeout, args.max_indegree, args.min_improvement)
        matlab_path = compare_root / "results" / "GTT" / f"{input_path.stem}_100_200_GTT_BIC.csv"
    else:
        result = bs_search(
            scorer,
            timeout,
            args.max_indegree,
            args.restarts,
            args.perturbations,
            args.seed,
            args.min_improvement,
        )
        matlab_path = compare_root / "results" / "BS" / f"{input_path.stem}_100_200_BayesianSearch.csv"

    names = [str(column) for column in frame.columns]
    edges = [(names[parent], names[child]) for parent, child in zip(*np.nonzero(result.dag))]
    output_engine = f"pgmpy_matlab_style_{algorithm}"
    files = write_outputs(args.output_dir, f"{input_path.stem}_trial{args.trial}", output_engine, names, edges, result.dag)

    return {
        "input": str(input_path),
        "dataset": input_path.stem,
        "algorithm": algorithm,
        "backend": "pgmpy",
        "search": "GTT" if algorithm == "gtt" else "BayesianSearch-style",
        "score": "bic-d",
        "trial": args.trial,
        "n_samples": int(frame.shape[0]),
        "n_nodes": int(frame.shape[1]),
        "max_indegree": args.max_indegree,
        "timeout": result.timeout,
        "elapsed_seconds": result.elapsed_seconds,
        "steps": result.steps,
        "edge_count": int(result.dag.sum()),
        "score_value": result.score,
        "metrics": strict_eval_against(result.dag, truth),
        "matlab": matlab_row(matlab_path, args.trial),
        "extra": result.extra,
        "files": files,
        **meta,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--algorithm", choices=["gtt", "bs", "both"], default="both")
    parser.add_argument("--trial", type=int, default=1)
    parser.add_argument("--mat-var", default=None)
    parser.add_argument("--network", default=None)
    parser.add_argument("--timeout-seconds", type=float, default=600)
    parser.add_argument("--max-indegree", type=int, default=None)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--perturbations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--min-improvement", type=float, default=1e-6)
    parser.add_argument("--csv-nodes-in-rows", action="store_true")
    parser.add_argument("--compare-root", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("results") / "open_source_bic")
    parser.add_argument("--no-cache", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    algorithms = ["gtt", "bs"] if args.algorithm == "both" else [args.algorithm]
    for input_path in args.inputs:
        for algorithm in algorithms:
            print(json.dumps(run_one(args, input_path, algorithm), ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
