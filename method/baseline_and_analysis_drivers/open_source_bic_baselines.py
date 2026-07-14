"""Run open-source BIC hill-climbing BN structure baselines.

The script supports the repository's MATLAB dataset format, where each .mat
file stores a 1 x trial cell array and each trial is a nodes x samples matrix.
It can also read ordinary CSV files with samples in rows.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.io import loadmat

from paper_style_bic_baselines import eval_against, matlab_row, normalize_network_name, parse_network


MISSING_VALUE = -9999


def load_frame(input_path: Path, mat_var: str | None, trial: int, csv_nodes_in_rows: bool) -> tuple[pd.DataFrame, dict]:
    if input_path.suffix.lower() == ".mat":
        mat = loadmat(input_path)
        keys = [key for key in mat if not key.startswith("__")]
        if not keys:
            raise ValueError(f"No MATLAB variables found in {input_path}")
        key = mat_var or keys[0]
        if key not in mat:
            raise ValueError(f"MAT variable {key!r} not found. Available variables: {keys}")

        value = mat[key]
        if getattr(value, "dtype", None) == object:
            flat = value.ravel()
            if trial < 1 or trial > len(flat):
                raise ValueError(f"Trial must be in [1, {len(flat)}], got {trial}")
            data = np.asarray(flat[trial - 1])
        else:
            data = np.asarray(value)

        if data.ndim != 2:
            raise ValueError(f"Expected a 2-D trial matrix, got shape {data.shape}")

        frame = pd.DataFrame(data.T, columns=[f"X{i + 1}" for i in range(data.shape[0])])
        meta = {
            "format": "mat",
            "mat_variable": key,
            "trial": trial,
            "mat_trial_shape": list(data.shape),
        }
    else:
        frame = pd.read_csv(input_path)
        if csv_nodes_in_rows:
            frame = frame.T.reset_index(drop=True).T
            frame.columns = [f"X{i + 1}" for i in range(frame.shape[1])]
        meta = {"format": "csv", "trial": None}

    frame = frame.replace(MISSING_VALUE, pd.NA).replace(str(MISSING_VALUE), pd.NA)
    for column in frame.columns:
        frame[column] = frame[column].astype("category")

    return frame, meta


def arcs_to_adjacency(edges: Iterable[tuple[str, str]], names: list[str]) -> np.ndarray:
    index = {name: i for i, name in enumerate(names)}
    adjacency = np.zeros((len(names), len(names)), dtype=int)
    for parent, child in edges:
        adjacency[index[str(parent)], index[str(child)]] = 1
    return adjacency


def run_pyagrum_bic(frame: pd.DataFrame, max_indegree: int | None, prior: str) -> tuple[list[tuple[str, str]], dict]:
    import pyagrum as gum

    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, newline="") as handle:
        temp_csv = Path(handle.name)
        frame.to_csv(handle, index=False)

    try:
        learner = gum.BNLearner(str(temp_csv))
        learner.useGreedyHillClimbing()
        learner.useScoreBIC()
        if prior == "smoothing":
            learner.useSmoothingPrior(1.0)
        elif prior == "bdeu":
            learner.useBDeuPrior(1.0)
        else:
            learner.useNoPrior()
        if max_indegree is not None:
            learner.setMaxIndegree(max_indegree)
        started = time.perf_counter()
        bn = learner.learnBN()
        elapsed = time.perf_counter() - started
        edges = sorted((bn.variable(tail).name(), bn.variable(head).name()) for tail, head in bn.arcs())
        return edges, {
            "library": "pyagrum",
            "version": gum.__version__,
            "search": "GreedyHillClimbing",
            "score": "BIC",
            "prior": prior,
            "elapsed_seconds": elapsed,
        }
    finally:
        temp_csv.unlink(missing_ok=True)


def run_pgmpy_bic(frame: pd.DataFrame, max_indegree: int | None) -> tuple[list[tuple[str, str]], dict]:
    import pgmpy

    started = time.perf_counter()
    try:
        from pgmpy.causal_discovery import HillClimbSearch

        estimator = HillClimbSearch(
            scoring_method="bic-d",
            max_indegree=max_indegree,
            return_type="dag",
            show_progress=False,
        )
        estimator.fit(frame)
        dag = estimator.causal_graph_
    except (ImportError, AttributeError):
        from pgmpy.estimators import HillClimbSearch

        search = HillClimbSearch(frame)
        dag = search.estimate(
            scoring_method="bic-d",
            max_indegree=max_indegree,
            show_progress=False,
        )

    elapsed = time.perf_counter() - started
    edges = sorted((str(parent), str(child)) for parent, child in dag.edges())
    return edges, {
        "library": "pgmpy",
        "version": pgmpy.__version__,
        "search": "HillClimbSearch",
        "score": "bic-d",
        "elapsed_seconds": elapsed,
    }


def write_outputs(output_dir: Path, stem: str, engine: str, names: list[str], edges: list[tuple[str, str]], adjacency: np.ndarray) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{stem}_{engine}_bic"
    edge_path = prefix.with_name(prefix.name + "_edges.csv")
    adjacency_path = prefix.with_name(prefix.name + "_adjacency.csv")
    adjacency_numeric_path = prefix.with_name(prefix.name + "_adjacency_numeric.csv")

    pd.DataFrame(edges, columns=["parent", "child"]).to_csv(edge_path, index=False)
    pd.DataFrame(adjacency, index=names, columns=names).to_csv(adjacency_path)
    pd.DataFrame(adjacency).to_csv(adjacency_numeric_path, index=False, header=False)

    return {
        "edges_csv": str(edge_path),
        "adjacency_csv": str(adjacency_path),
        "adjacency_numeric_csv": str(adjacency_numeric_path),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input .mat or .csv dataset.")
    parser.add_argument("--mat-var", default=None, help="MATLAB variable name. Defaults to the first non-metadata variable.")
    parser.add_argument("--trial", type=int, default=1, help="1-based trial index for .mat cell-array datasets.")
    parser.add_argument("--engine", choices=["pyagrum", "pgmpy", "both"], default="both")
    parser.add_argument("--max-indegree", type=int, default=None, help="Optional maximum number of parents per node.")
    parser.add_argument("--csv-nodes-in-rows", action="store_true", help="Use when CSV rows are nodes and columns are samples.")
    parser.add_argument("--output-dir", type=Path, default=Path("results") / "open_source_bic")
    parser.add_argument("--compare-root", type=Path, default=Path("."), help="Repository root for true networks and MATLAB result CSVs.")
    parser.add_argument("--network", default=None, help="Benchmark network name for structural evaluation. Defaults to the MAT variable/input stem.")
    parser.add_argument("--pyagrum-prior", choices=["none", "smoothing", "bdeu"], default="none", help="Optional pyAgrum prior. Default preserves no-prior BIC; smoothing/BDeu can avoid sparse-count failures.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    frame, meta = load_frame(args.input, args.mat_var, args.trial, args.csv_nodes_in_rows)
    names = [str(column) for column in frame.columns]
    stem = args.input.stem if meta["format"] == "csv" else f"{args.input.stem}_trial{args.trial}"

    engines = ["pyagrum", "pgmpy"] if args.engine == "both" else [args.engine]
    runners = {
        "pyagrum": run_pyagrum_bic,
        "pgmpy": run_pgmpy_bic,
    }

    summary = {
        "input": str(args.input),
        "n_samples": int(frame.shape[0]),
        "n_nodes": int(frame.shape[1]),
        "columns": names,
        "scoring": "BIC",
        "max_indegree": args.max_indegree,
        "pyagrum_prior": args.pyagrum_prior,
        **meta,
        "results": {},
    }

    truth = None
    matlab_comparison = {}
    if meta["format"] == "mat":
        compare_root = args.compare_root.resolve()
        network_name = args.network or normalize_network_name(args.mat_var or args.input.stem)
        data_for_sizes = np.asarray(frame.astype(int).to_numpy().T, dtype=np.int64)
        truth, _node_sizes = parse_network(network_name, compare_root, data_for_sizes)
        matlab_comparison = {
            "gtt": matlab_row(compare_root / "results" / "GTT" / f"{args.input.stem}_100_200_GTT_BIC.csv", args.trial),
            "bs": matlab_row(compare_root / "results" / "BS" / f"{args.input.stem}_100_200_BayesianSearch.csv", args.trial),
        }
        summary["comparison"] = {
            "network": network_name,
            "matlab": matlab_comparison,
        }

    failures = 0
    for engine in engines:
        try:
            if engine == "pyagrum":
                edges, engine_meta = run_pyagrum_bic(frame, args.max_indegree, args.pyagrum_prior)
            else:
                edges, engine_meta = runners[engine](frame, args.max_indegree)
            adjacency = arcs_to_adjacency(edges, names)
            output_engine = engine if engine != "pyagrum" or args.pyagrum_prior == "none" else f"{engine}_{args.pyagrum_prior}"
            files = write_outputs(args.output_dir, stem, output_engine, names, edges, adjacency)
            summary["results"][engine] = {
                **engine_meta,
                "edge_count": len(edges),
                "metrics": eval_against(adjacency, truth) if truth is not None else {},
                "edges": edges,
                "files": files,
            }
        except Exception as exc:  # Keep the other engine usable if one package has stricter data requirements.
            failures += 1
            summary["results"][engine] = {"error": f"{type(exc).__name__}: {exc}"}

    print(json.dumps(summary, indent=2))
    return 1 if failures == len(engines) else 0


if __name__ == "__main__":
    sys.exit(main())
