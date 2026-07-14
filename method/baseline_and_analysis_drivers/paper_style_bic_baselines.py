"""Paper-style GTT and Bayesian Search baselines with BIC scoring.

This script is intended for verification against the MATLAB reproductions in
results/GTT and results/BS. It mirrors the project implementations:

- GTT_BIC: empty graph -> greedy thickening by BIC -> greedy thinning by BIC.
- BayesianSearch: maximum-branching initialization, greedy add/delete/reverse
  local search, then random perturbation + local-search restarts.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat


MISSING_VALUE = -9999
NETWORK_FILE = {
    "asia": "mk_asia_bnet.m",
    "alarm": "mk_alarm_bnet.m",
    "ins": "mk_insur_bnet.m",
    "insurance": "mk_insur_bnet.m",
    "water": "mk_water_bnet.m",
    "hailfinder": "mk_hailfinder_bnet.m",
    "hepar": "mk_hepar2_bnet.m",
    "hepar ii": "mk_hepar2_bnet.m",
    "win95pts": "mk_win95pts_bnet.m",
    "and": "mk_andes_bnet.m",
    "andes": "mk_andes_bnet.m",
}


@dataclass
class SearchResult:
    dag: np.ndarray
    score: float
    steps: int
    timeout: bool
    elapsed_seconds: float
    extra: dict


class Timeout:
    def __init__(self, seconds: float | None):
        self.seconds = seconds
        self.started = time.perf_counter()

    def expired(self) -> bool:
        return self.seconds is not None and (time.perf_counter() - self.started) >= self.seconds


class BICScorer:
    def __init__(self, data_nodes_by_samples: np.ndarray, node_sizes: np.ndarray):
        self.data = data_nodes_by_samples.astype(np.int64, copy=False)
        self.ns = node_sizes.astype(np.int64, copy=False)
        self.n_nodes = self.data.shape[0]
        self.cache: dict[tuple[int, tuple[int, ...]], float] = {}

    def local_score(self, child: int, parents: list[int] | tuple[int, ...] | np.ndarray) -> float:
        parent_tuple = tuple(sorted(int(p) for p in parents))
        key = (int(child), parent_tuple)
        if key in self.cache:
            return self.cache[key]

        family = (child,) + parent_tuple
        family_data = self.data[list(family), :]
        available = np.all(family_data != MISSING_VALUE, axis=0)
        child_values = self.data[child, available] - 1
        sample_count = int(child_values.size)
        r = int(self.ns[child])
        q = int(np.prod(self.ns[list(parent_tuple)])) if parent_tuple else 1

        if sample_count == 0:
            score = -math.inf
            self.cache[key] = score
            return score

        if parent_tuple:
            parent_data = self.data[list(parent_tuple), :][:, available] - 1
            multipliers = np.concatenate(([1], np.cumprod(self.ns[list(parent_tuple[:-1])]))) if len(parent_tuple) > 1 else np.array([1])
            parent_index = (parent_data.T * multipliers).sum(axis=1).astype(np.int64)
        else:
            parent_index = np.zeros(sample_count, dtype=np.int64)

        valid = (child_values >= 0) & (child_values < r) & (parent_index >= 0) & (parent_index < q)
        flat = parent_index[valid] * r + child_values[valid]
        counts = np.bincount(flat, minlength=q * r).reshape(q, r).astype(float)

        # Match the packaged MATLAB baseline: score_family(..., 'bic') builds a
        # tabular_CPD with prior_type='dirichlet', dirichlet_weight=1. BNT then
        # learns a smoothed CPT and evaluates the data log-likelihood under that
        # CPT before applying the BIC penalty.
        alpha = 1.0 / (q * r)
        smoothed = counts + alpha
        cpt = smoothed / smoothed.sum(axis=1, keepdims=True)
        nonzero = counts > 0
        log_likelihood = float(np.sum(counts[nonzero] * np.log(cpt[nonzero])))
        nparams = (r - 1) * q
        score = log_likelihood - 0.5 * nparams * math.log(sample_count)
        self.cache[key] = score
        return score

    def dag_score(self, dag: np.ndarray) -> float:
        total = 0.0
        for child in range(self.n_nodes):
            total += self.local_score(child, np.flatnonzero(dag[:, child]))
        return total


def load_trial(mat_path: Path, trial: int) -> tuple[np.ndarray, str]:
    mat = loadmat(mat_path)
    keys = [key for key in mat if not key.startswith("__")]
    if not keys:
        raise ValueError(f"No MATLAB variables found in {mat_path}")
    key = keys[0]
    value = mat[key]
    if getattr(value, "dtype", None) == object:
        flat = value.ravel()
        if trial < 1 or trial > len(flat):
            raise ValueError(f"Trial must be in [1, {len(flat)}], got {trial}")
        data = np.asarray(flat[trial - 1], dtype=np.int64)
    else:
        data = np.asarray(value, dtype=np.int64)
    if data.ndim != 2:
        raise ValueError(f"Expected nodes x samples matrix, got {data.shape}")
    return data, key


def normalize_network_name(name: str) -> str:
    stem = re.sub(r"\d+$", "", name).lower()
    if stem == "ins":
        return "insurance"
    if stem == "hepar":
        return "hepar"
    return stem


def parse_network(network: str, root: Path, fallback_data: np.ndarray | None = None) -> tuple[np.ndarray | None, np.ndarray]:
    key = normalize_network_name(network)
    filename = NETWORK_FILE.get(key)
    if filename is None:
        if fallback_data is None:
            raise ValueError(f"Unknown network {network!r}")
        return None, fallback_node_sizes(fallback_data)
    text = (root / "BS_GTT_standalone_package_20260618" / "networks" / filename).read_text(encoding="utf-8", errors="ignore")
    env = parse_scalar_assignments(text)
    env.update(parse_node_struct(text))
    n = parse_node_count(text, env, fallback_data)
    adjacency = np.zeros((n, n), dtype=np.int8)

    for var_name in ("dag", "adjacency"):
        pattern = re.compile(rf"{var_name}\((.+?),\s*(.+?)\)\s*=\s*1\s*;", re.S)
        for left, right in pattern.findall(text):
            try:
                parents = parse_index_expr(left, env)
                children = parse_index_expr(right, env)
            except ValueError:
                continue
            for parent in parents:
                for child in children:
                    if 1 <= parent <= n and 1 <= child <= n:
                        adjacency[parent - 1, child - 1] = 1

    sizes = parse_node_sizes(text, n)
    if sizes is None:
        sizes = fallback_node_sizes(fallback_data) if fallback_data is not None else np.ones(n, dtype=np.int64) * 2
    return adjacency, sizes


def parse_scalar_assignments(text: str) -> dict[str, int]:
    env: dict[str, int] = {}
    for name, value in re.findall(r"^\s*([A-Za-z]\w*)\s*=\s*(\d+)\s*;", text, flags=re.M):
        env[name] = int(value)
    return env


def parse_node_struct(text: str) -> dict[str, int]:
    env: dict[str, int] = {}
    match = re.search(r"node\s*=\s*struct\((.*?)\)\s*;", text, flags=re.S)
    if not match:
        return env
    for name, value in re.findall(r"'([^']+)'\s*,\s*(\d+)", match.group(1)):
        env["node." + name] = int(value)
    return env


def parse_node_count(text: str, env: dict[str, int], fallback_data: np.ndarray | None) -> int:
    for pattern in (r"\bN\s*=\s*(\d+)\s*;", r"\bn\s*=\s*(\d+)\s*;", r"(?:dag|adjacency)\s*=\s*zeros\((\d+)"):
        match = re.search(pattern, text)
        if match:
            return int(match.group(1))
    if env:
        return max(env.values())
    if fallback_data is not None:
        return int(fallback_data.shape[0])
    raise ValueError("Cannot determine node count")


def parse_index_expr(expr: str, env: dict[str, int]) -> list[int]:
    expr = expr.strip()
    if expr.startswith("[") and expr.endswith("]"):
        expr = expr[1:-1]
    expr = expr.replace(",", " ")
    values: list[int] = []
    for token in expr.split():
        token = token.strip()
        if not token:
            continue
        if token.isdigit():
            values.append(int(token))
        elif token in env:
            values.append(env[token])
        else:
            raise ValueError(token)
    return values


def parse_node_sizes(text: str, n: int) -> np.ndarray | None:
    match = re.search(r"node_sizes\s*=\s*(\d+)\s*\*\s*ones\(1\s*,\s*N\)\s*;", text)
    if match:
        sizes = np.ones(n, dtype=np.int64) * int(match.group(1))
        for idx, value in re.findall(r"node_sizes\((\d+)\)\s*=\s*(\d+)\s*;", text):
            sizes[int(idx) - 1] = int(value)
        return sizes

    match = re.search(r"mk_bnet\([^,]+,\s*\[([^\]]+)\]", text, flags=re.S)
    if match:
        numbers = [int(x) for x in re.findall(r"\d+", match.group(1))]
        if len(numbers) == n:
            return np.asarray(numbers, dtype=np.int64)

    if re.search(r"ns\s*=\s*arity\s*\*\s*ones", text):
        return np.ones(n, dtype=np.int64) * 2
    return None


def fallback_node_sizes(data: np.ndarray | None) -> np.ndarray:
    if data is None:
        raise ValueError("Cannot infer node sizes without data")
    masked = np.where(data == MISSING_VALUE, 0, data)
    return masked.max(axis=1).astype(np.int64)


def has_directed_path(dag: np.ndarray, source: int, target: int) -> bool:
    if source == target:
        return True
    visited = np.zeros(dag.shape[0], dtype=bool)
    stack = [int(source)]
    while stack:
        node = stack.pop()
        if visited[node]:
            continue
        visited[node] = True
        children = np.flatnonzero(dag[node, :])
        if np.any(children == target):
            return True
        stack.extend(int(child) for child in children if not visited[child])
    return False


def gtt_search(
    scorer: BICScorer,
    timeout: Timeout,
    max_indegree: int | None = None,
    min_improvement: float = 0.0,
) -> SearchResult:
    start = time.perf_counter()
    n = scorer.n_nodes
    dag = np.zeros((n, n), dtype=np.int8)
    local = np.array([scorer.local_score(child, []) for child in range(n)])
    score = float(local.sum())
    thick_steps = thin_steps = 0

    while not timeout.expired():
        best = None
        best_delta = min_improvement
        for parent in range(n):
            for child in range(n):
                if parent == child or dag[parent, child] or dag[child, parent]:
                    continue
                old_parents = list(np.flatnonzero(dag[:, child]))
                if has_directed_path(dag, child, parent):
                    continue
                new_score = scorer.local_score(child, old_parents + [parent])
                delta = new_score - local[child]
                if delta > best_delta:
                    best_delta = delta
                    best = (parent, child, new_score)
        if best is None:
            break
        parent, child, new_score = best
        dag[parent, child] = 1
        local[child] = new_score
        score += best_delta
        thick_steps += 1

    while not timeout.expired():
        best = None
        best_delta = min_improvement
        for parent, child in zip(*np.nonzero(dag)):
            old_parents = list(np.flatnonzero(dag[:, child]))
            new_parents = [p for p in old_parents if p != parent]
            new_score = scorer.local_score(child, new_parents)
            delta = new_score - local[child]
            if delta > best_delta:
                best_delta = delta
                best = (int(parent), int(child), new_score)
        if best is None:
            break
        parent, child, new_score = best
        dag[parent, child] = 0
        local[child] = new_score
        score += best_delta
        thin_steps += 1

    dag, score, local, forced_steps = enforce_parent_limit_deletions(dag, score, local, scorer, timeout, max_indegree)
    thin_steps += forced_steps

    return SearchResult(dag, score, thick_steps + thin_steps, timeout.expired(), time.perf_counter() - start, {"thick_steps": thick_steps, "thin_steps": thin_steps})


def bs_search(
    scorer: BICScorer,
    timeout: Timeout,
    max_indegree: int | None = None,
    restarts: int = 2,
    perturbations: int = 3,
    seed: int = 1,
    min_improvement: float = 0.0,
) -> SearchResult:
    start = time.perf_counter()
    rng = np.random.default_rng(seed)
    dag, local = maximum_branching_start(scorer)
    score = float(local.sum())
    dag, score, local, steps = hill_climb(dag, score, local, scorer, timeout, max_indegree, min_improvement)
    dag, score, local, limit_steps = enforce_parent_limit_deletions(dag, score, local, scorer, timeout, max_indegree)
    best_dag = dag.copy()
    best_score = score
    total_steps = steps + limit_steps
    perturb_steps = 0
    completed = 0

    for _ in range(restarts):
        if timeout.expired():
            break
        dag, score, local, used = perturb_dag(dag, score, local, scorer, timeout, rng, max_indegree, perturbations)
        perturb_steps += used
        dag, score, local, steps = hill_climb(dag, score, local, scorer, timeout, max_indegree, min_improvement)
        dag, score, local, limit_steps = enforce_parent_limit_deletions(dag, score, local, scorer, timeout, max_indegree)
        total_steps += steps + limit_steps
        completed += 1
        if score > best_score:
            best_score = score
            best_dag = dag.copy()

    return SearchResult(best_dag, best_score, total_steps, timeout.expired(), time.perf_counter() - start, {"perturbation_steps": perturb_steps, "completed_restarts": completed})


def maximum_branching_start(scorer: BICScorer) -> tuple[np.ndarray, np.ndarray]:
    n = scorer.n_nodes
    empty = np.array([scorer.local_score(child, []) for child in range(n)])
    weights = np.full((n + 1, n + 1), -math.inf)
    root = n
    weights[root, :n] = 0.0
    single = np.zeros((n, n))
    for parent in range(n):
        for child in range(n):
            if parent == child:
                continue
            single[parent, child] = scorer.local_score(child, [parent])
            weights[parent, child] = single[parent, child] - empty[child]
    parents = maximum_branching_parents(weights, root)
    dag = np.zeros((n, n), dtype=np.int8)
    local = empty.copy()
    for child in range(n):
        parent = parents[child]
        if parent >= 0 and parent != root:
            dag[parent, child] = 1
            local[child] = single[parent, child]
    return dag, local


def maximum_branching_parents(weights: np.ndarray, root: int) -> np.ndarray:
    n = weights.shape[0]
    parent = np.full(n, -1, dtype=int)
    for child in range(n):
        if child == root:
            continue
        parent[child] = int(np.argmax(weights[:, child]))
    cycle = find_parent_cycle(parent, root)
    if not cycle:
        return parent
    in_cycle = np.zeros(n, dtype=bool)
    in_cycle[cycle] = True
    outside = [i for i in range(n) if not in_cycle[i]]
    super_node = len(outside)
    old_to_new = np.full(n, -1, dtype=int)
    for idx, old in enumerate(outside):
        old_to_new[old] = idx
    for old in cycle:
        old_to_new[old] = super_node
    new_to_old = outside + [-1]
    new_root = old_to_new[root]
    new_weights = np.full((super_node + 1, super_node + 1), -math.inf)
    enter_target = np.full(super_node + 1, -1, dtype=int)
    exit_source = np.full(super_node + 1, -1, dtype=int)
    for u in outside:
        for v in outside:
            new_weights[old_to_new[u], old_to_new[v]] = weights[u, v]
    for u in outside:
        vals = [(weights[u, v] - weights[parent[v], v], v) for v in cycle]
        best_weight, best_v = max(vals, key=lambda x: x[0])
        new_weights[old_to_new[u], super_node] = best_weight
        enter_target[old_to_new[u]] = best_v
    for v in outside:
        vals = [(weights[u, v], u) for u in cycle]
        best_weight, best_u = max(vals, key=lambda x: x[0])
        new_weights[super_node, old_to_new[v]] = best_weight
        exit_source[old_to_new[v]] = best_u
    contracted_parent = maximum_branching_parents(new_weights, new_root)
    expanded = parent.copy()
    expanded[root] = -1
    for old in outside:
        np_child = old_to_new[old]
        np_parent = contracted_parent[np_child]
        if np_parent == super_node:
            expanded[old] = exit_source[np_child]
        elif np_parent < 0:
            expanded[old] = -1
        else:
            expanded[old] = new_to_old[np_parent]
    super_parent = contracted_parent[super_node]
    if super_parent >= 0:
        old_parent = new_to_old[super_parent]
        replaced_child = enter_target[super_parent]
        expanded[replaced_child] = old_parent
    return expanded


def find_parent_cycle(parent: np.ndarray, root: int) -> list[int]:
    for start in range(len(parent)):
        if start == root:
            continue
        order: dict[int, int] = {}
        path: list[int] = []
        node = start
        while node >= 0 and node != root and node not in order:
            order[node] = len(path)
            path.append(node)
            node = int(parent[node])
        if node >= 0 and node != root and node in order:
            return path[order[node] :]
    return []


def hill_climb(
    dag: np.ndarray,
    score: float,
    local: np.ndarray,
    scorer: BICScorer,
    timeout: Timeout,
    max_indegree: int | None,
    min_improvement: float = 0.0,
) -> tuple[np.ndarray, float, np.ndarray, int]:
    steps = 0
    while not timeout.expired():
        move = best_hill_move(dag, local, scorer, timeout, max_indegree, min_improvement)
        if move is None:
            break
        dag, score, local = apply_move(dag, score, local, move)
        steps += 1
    return dag, score, local, steps


def best_hill_move(
    dag: np.ndarray,
    local: np.ndarray,
    scorer: BICScorer,
    timeout: Timeout,
    max_indegree: int | None,
    min_improvement: float = 0.0,
):
    n = dag.shape[0]
    best = None
    best_delta = min_improvement
    for parent in range(n):
        for child in range(n):
            if timeout.expired() or parent == child:
                continue
            if not dag[parent, child] and not dag[child, parent]:
                parents = list(np.flatnonzero(dag[:, child]))
                if not has_directed_path(dag, child, parent):
                    new_child_score = scorer.local_score(child, parents + [parent])
                    delta = new_child_score - local[child]
                    if delta > best_delta:
                        best_delta = delta
                        best = ("add", parent, child, delta, new_child_score, None)
            elif dag[parent, child]:
                parents = list(np.flatnonzero(dag[:, child]))
                new_parents = [p for p in parents if p != parent]
                new_child_score = scorer.local_score(child, new_parents)
                delta = new_child_score - local[child]
                if delta > best_delta:
                    best_delta = delta
                    best = ("delete", parent, child, delta, new_child_score, None)
                parent_parents = list(np.flatnonzero(dag[:, parent]))
                candidate = dag.copy()
                candidate[parent, child] = 0
                if not has_directed_path(candidate, parent, child):
                    new_parent_score = scorer.local_score(parent, parent_parents + [child])
                    delta = new_child_score + new_parent_score - local[child] - local[parent]
                    if delta > best_delta:
                        best_delta = delta
                        best = ("reverse", parent, child, delta, new_child_score, new_parent_score)
    return best


def enforce_parent_limit_deletions(
    dag: np.ndarray,
    score: float,
    local: np.ndarray,
    scorer: BICScorer,
    timeout: Timeout,
    max_indegree: int | None,
) -> tuple[np.ndarray, float, np.ndarray, int]:
    steps = 0
    while not timeout.expired():
        move = forced_parent_limit_delete(dag, local, scorer, max_indegree)
        if move is None:
            break
        dag, score, local = apply_move(dag, score, local, move)
        steps += 1
    return dag, score, local, steps


def forced_parent_limit_delete(dag: np.ndarray, local: np.ndarray, scorer: BICScorer, max_indegree: int | None):
    if max_indegree is None:
        return None
    best = None
    best_delta = -math.inf
    for child in range(dag.shape[1]):
        parents = list(np.flatnonzero(dag[:, child]))
        if len(parents) <= max_indegree:
            continue
        for parent in parents:
            new_parents = [p for p in parents if p != parent]
            new_child_score = scorer.local_score(child, new_parents)
            delta = new_child_score - local[child]
            if delta > best_delta:
                best_delta = delta
                best = ("delete", int(parent), int(child), delta, new_child_score, None)
    return best


def perturb_dag(dag: np.ndarray, score: float, local: np.ndarray, scorer: BICScorer, timeout: Timeout, rng: np.random.Generator, max_indegree: int | None, perturbations: int):
    used = 0
    n = dag.shape[0]
    for _ in range(perturbations):
        if timeout.expired():
            break
        move = None
        for _attempt in range(max(1000, 20 * n * n)):
            parent = int(rng.integers(n))
            child = int(rng.integers(n))
            if parent == child:
                continue
            move_type = int(rng.integers(3))
            if move_type == 0 and not dag[parent, child] and not dag[child, parent]:
                parents = list(np.flatnonzero(dag[:, child]))
                if not has_directed_path(dag, child, parent):
                    new_child_score = scorer.local_score(child, parents + [parent])
                    move = ("add", parent, child, new_child_score - local[child], new_child_score, None)
                    break
            elif move_type == 1 and dag[parent, child]:
                parents = list(np.flatnonzero(dag[:, child]))
                new_child_score = scorer.local_score(child, [p for p in parents if p != parent])
                move = ("delete", parent, child, new_child_score - local[child], new_child_score, None)
                break
            elif move_type == 2 and dag[parent, child]:
                parent_parents = list(np.flatnonzero(dag[:, parent]))
                candidate = dag.copy()
                candidate[parent, child] = 0
                if not has_directed_path(candidate, parent, child):
                    child_parents = [p for p in np.flatnonzero(dag[:, child]) if p != parent]
                    new_child_score = scorer.local_score(child, child_parents)
                    new_parent_score = scorer.local_score(parent, parent_parents + [child])
                    move = ("reverse", parent, child, new_child_score + new_parent_score - local[child] - local[parent], new_child_score, new_parent_score)
                    break
        if move is None:
            break
        dag, score, local = apply_move(dag, score, local, move)
        used += 1
    return dag, score, local, used


def apply_move(dag: np.ndarray, score: float, local: np.ndarray, move):
    kind, parent, child, delta, new_child_score, new_parent_score = move
    dag = dag.copy()
    local = local.copy()
    if kind == "add":
        dag[parent, child] = 1
        local[child] = new_child_score
    elif kind == "delete":
        dag[parent, child] = 0
        local[child] = new_child_score
    elif kind == "reverse":
        dag[parent, child] = 0
        dag[child, parent] = 1
        local[child] = new_child_score
        local[parent] = new_parent_score
    return dag, score + delta, local


def eval_against(pred: np.ndarray, truth: np.ndarray | None) -> dict:
    if truth is None:
        return {}
    n = truth.shape[0]
    tp = fn = fp = tn = 0
    for i in range(n):
        for j in range(i, n):
            t_ij, t_ji = truth[i, j], truth[j, i]
            p_ij, p_ji = pred[i, j], pred[j, i]
            if t_ij or t_ji:
                if (t_ij and p_ij) or (t_ji and p_ji):
                    tp += 1
                elif (t_ij and p_ji) or (t_ji and p_ij):
                    fp += 1
                else:
                    fn += 1
            else:
                if p_ij or p_ji:
                    fp += 1
                else:
                    tn += 1
    sensitivity = tp / (tp + fn) if (tp + fn) else 1.0
    precision = tp / (tp + fp) if (tp + fp) else 1.0
    f1 = 2 * sensitivity * precision / (sensitivity + precision) if (sensitivity + precision) else 0.0
    specificity = tn / (tn + fp) if (tn + fp) else 0.0
    shd = fn + fp
    return {
        "f1": f1,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "shd": shd,
        "tp": tp,
        "fn": fn,
        "fp": fp,
        "tn": tn,
    }


def matlab_row(path: Path, trial: int) -> dict | None:
    if not path.exists():
        return None
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    for row in rows:
        if row and row[0].strip().isdigit() and int(row[0]) == trial:
            vals = [x.strip() for x in row]
            if len(vals) >= 9 and vals[1].lower() not in {"timeout", "memoryoverflow"}:
                return {
                    "f1": float(vals[1]),
                    "sensitivity": float(vals[2]),
                    "specificity": float(vals[3]),
                    "precision": float(vals[4]),
                    "shd": float(vals[5]),
                    "score": float(vals[6]),
                    "elapsed_seconds": float(vals[7]),
                    "steps": float(vals[8]),
                }
            return {"status": vals[1], "elapsed_seconds": float(vals[2]) if len(vals) > 2 else None}
    return None


def run_one(args: argparse.Namespace, mat_path: Path, algorithm: str) -> dict:
    data, mat_var = load_trial(mat_path, args.trial)
    network = args.network or normalize_network_name(mat_var)
    truth, node_sizes = parse_network(network, args.root, data)
    scorer = BICScorer(data, node_sizes)
    timeout = Timeout(args.timeout_seconds)
    if algorithm == "gtt":
        result = gtt_search(scorer, timeout, args.max_indegree)
        matlab_path = args.root / "results" / "GTT" / f"{mat_path.stem}_100_200_GTT_BIC.csv"
    else:
        result = bs_search(scorer, timeout, args.max_indegree, args.restarts, args.perturbations, args.seed)
        matlab_path = args.root / "results" / "BS" / f"{mat_path.stem}_100_200_BayesianSearch.csv"
    metrics = eval_against(result.dag, truth)
    row = {
        "dataset": mat_path.stem,
        "algorithm": algorithm,
        "trial": args.trial,
        "n_nodes": int(data.shape[0]),
        "n_samples": int(data.shape[1]),
        "score": result.score,
        "elapsed_seconds": result.elapsed_seconds,
        "steps": result.steps,
        "edge_count": int(result.dag.sum()),
        "timeout": result.timeout,
        **metrics,
        **result.extra,
        "matlab": matlab_row(matlab_path, args.trial),
    }
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--algorithm", choices=["gtt", "bs", "both"], default="both")
    parser.add_argument("--trial", type=int, default=1)
    parser.add_argument("--network", default=None)
    parser.add_argument("--timeout-seconds", type=float, default=600)
    parser.add_argument("--max-indegree", type=int, default=None)
    parser.add_argument("--restarts", type=int, default=2)
    parser.add_argument("--perturbations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.root = args.root.resolve()
    algorithms = ["gtt", "bs"] if args.algorithm == "both" else [args.algorithm]
    results = []
    for input_path in args.inputs:
        for algorithm in algorithms:
            results.append(run_one(args, input_path, algorithm))
            print(json.dumps(results[-1], ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
