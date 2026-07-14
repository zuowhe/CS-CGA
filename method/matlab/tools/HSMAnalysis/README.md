# HSM Cache Behavior and Tradeoff Analysis

This folder contains the experiment code added to address reviewer comment 5:

> The evidence for HSM is incomplete. The paper emphasizes efficiency gains
> related to memoization, yet the experiments mainly report runtime
> improvements and do not provide enough direct analysis of cache hit rate,
> memory consumption, or the tradeoff between cache size and performance.

The goal of this experiment is not to repeat the full performance comparison.
Instead, it directly analyzes the cache behavior of HSM and its time-memory
tradeoff.

## Experiment Scope

Each benchmark network is evaluated at the 1000-sample setting and is run once.
This is intentional: the purpose is cache behavior analysis, not a new
statistical performance comparison.

The default datasets in `run_HSM_runtime_ablation.m` are:

- `Asia-1000`
- `INS-1000`
- `Water-1000`
- `Alarm-1000`
- `Hailfinder-1000`
- `HEPAR-1000`
- `Win95pts-1000`
- `AND-1000`

The run uses a fixed base seed so that cache configurations are compared under
the same search setting as much as possible.

## Cache Variants

The variants are defined in `HSM_analysis_default_variants.m`.

`BNT-Cache`:
A standard local score cache based on the BNT scoring cache. It caches local
family scores by child node and parent set.

`GB-Cache`:
A hybrid baseline that combines a global DAG topology hash with the standard
BNT local score cache. The global hash avoids rescoring a complete DAG when the
same topology appears again.

`Large-Hash`:
A large hash-table cache mode. It stores both global DAG scores and local
parent-set scores in hash maps and uses a 4 GB target capacity marker. It is
included as a contrast to show whether a coarse large-cache strategy is more
memory efficient than HSM. The code records the actual number of entries and an
estimated memory footprint rather than preallocating the full 4 GB.

`HSM-50n`, `HSM-100n`, `HSM-200n`, `HSM-500n`:
HSM variants with different per-node parent-set cache capacities. Here `n` is
the number of nodes in the network. For example, `HSM-200n` sets the maximum
parent-set entries for each node to:

```matlab
max_cache_size = 200 * BN_NodesNum
```

`HSM-200n` is the default setting used by the proposed method.

`HSM-large`:
An approximately unbounded parent-set cache mode using a very large per-node
capacity. It is used to check whether increasing the HSM cache beyond the
default setting still provides meaningful runtime improvement.

## Recorded Metrics

The output CSV records both performance and cache-behavior metrics:

- `runtimeSec`
- `f1`
- `shd`
- `globalQueries`
- `globalHits`
- `globalHitRate`
- `localQueries`
- `localHits`
- `localHitRate`
- `globalEntries`
- `localEntries`
- `overwrites`
- `bntCacheCapacity`
- `parentCacheSize`
- `largeHashTargetMB`
- `estimatedCacheMB`
- `matlabUsedMB`

The key metrics for responding to the reviewer are:

```text
globalHitRate = globalHits / globalQueries
localHitRate  = localHits / localQueries
```

The memory numbers are estimates. They are intended to represent the cache
footprint used by the different strategies, not an exact MATLAB object memory
profile.

For HSM, the local cache estimate is based on the number of cached parent sets:

```text
local cache memory ~= cached parent sets * (node mask + score value)
```

For global topology caches, the estimate uses stored DAG-key characters plus
stored scores:

```text
global cache memory ~= DAG key characters + cached score values
```

The code also reports `matlabUsedMB` from MATLAB's `memory` function when that
information is available.

## Expected Analysis Logic

This experiment supports two related claims.

First, the structural comparison (`BNT-Cache`, `GB-Cache`, `Large-Hash`, and
default `HSM-200n`) shows whether HSM is more practical than simply using a
large hash table. If `Large-Hash` requires much more estimated memory while
providing limited additional runtime improvement, then HSM's hierarchical
design is better justified.

Second, the HSM cache-size sensitivity analysis (`HSM-50n`, `HSM-100n`,
`HSM-200n`, `HSM-500n`, and `HSM-large`) shows the tradeoff between cache
capacity, hit rate, runtime, and memory usage.

A useful conclusion to check in the results is:

> As the parent-set cache capacity increases, the local hit rate improves and
> runtime decreases, but the marginal gain becomes small after `HSM-200n`.
> The large-cache setting further reduces overwrites, but provides limited
> additional runtime improvement while requiring more memory. Therefore,
> `HSM-200n` offers a practical tradeoff between acceleration and memory
> consumption.

F1 and SHD are recorded mainly to verify that cache settings do not materially
change the learned structures. HSM changes score reuse, not the search
operators, so F1 and SHD should remain unchanged or nearly unchanged under a
fixed seed.

## How to Run

From MATLAB:

```matlab
cd('path/to/CS_CGA_open_source/src/matlab')
setup_CS_CGA_paths
run_HSM_runtime_ablation
```

From a Windows command prompt:

```bat
cd /d "path\to\CS_CGA_open_source\src\matlab"
matlab -batch "run_HSM_runtime_ablation"
```

The summary file is written to:

```text
results/HSM_analysis_YYYYMMDD_HHMMSS/hsm_analysis_summary.csv
```

## Code Map

- `run_HSM_runtime_ablation.m`: top-level experiment script.
- `HSM_analysis_default_variants.m`: cache variant definitions.
- `HSM_analysis_init_state.m`: initializes cache state and capacities.
- `HSM_analysis_score_population.m`: scores populations and records cache hits.
- `HSM_analysis_cache_stats.m`: converts cache state into report metrics.
- `HSM_analysis_process.m`: CS-CGA run with instrumented scoring cache.
- `setup_HSM_analysis_paths.m`: adds the required paths for this analysis.
