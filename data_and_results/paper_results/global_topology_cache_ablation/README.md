# Global Topology Cache Ablation Results

These files support the individual-level global topology cache analysis in
Section 5.5.2 of the manuscript. Each configuration is evaluated on the eight
benchmark BNs with 1000 samples and 200 generations.

## Configurations

- `HSM`: the complete proposed strategy, including complete-DAG score reuse
  through hash-based global topology lookup and local parent-set caching.
- `NoReuse`: retains local parent-set caching and historical DAG keys for
  duplicate handling, but does not reuse cached complete-DAG scores.
- `NoHash`: retains the same complete-DAG score cache as HSM but retrieves
  scores through sequential DAG matching rather than hash-based lookup.

## Files

- `global_topology_cache_ablation_table.csv`: the eight-row summary used in
  the manuscript table. The cache-hit and cache-size fields describe HSM.
- `hsm_global_topology_cache.csv`: selected `HSM-paper` rows from the HSM
  analysis output.
- `noreuse_global_topology_cache.csv`: NoReuse results.
- `nohash_global_topology_cache.csv`: NoHash results.

`hsmGlobalCacheKB` estimates the storage of cached DAG keys and complete-DAG
scores. It excludes the local parent-set cache and is not MATLAB process memory.
