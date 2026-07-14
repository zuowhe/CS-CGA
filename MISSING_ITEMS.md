# Publication Readiness Checklist

## Must be completed before publication

1. Add a top-level `LICENSE` after confirming the intended license for the CS-CGA
   code. Do not redistribute third-party packages unless their licenses permit it.
2. Add a public-facing README with exact commands for: CS-CGA experiments,
   baseline runs, result aggregation, Wilcoxon tests, and every figure.
3. Record the MATLAB release, BNT version and installation procedure, Python
   version, and exact `pgmpy`, `numpy`, `pandas`, `matplotlib`, and `seaborn`
   versions. A basic `plotting/requirements.txt` is included, but the versions
   must be pinned after a clean-environment reproduction test.
4. Make the latest HSM rerun scripts self-contained. The files in
   `method/hsm_analysis_scripts/` currently expect `tools/HSMAnalysis` and
   `tools/AESLHSMAnalysis` next to the scripts. Their required latest support
   code has not yet been collected into this staging package. The bundled
   `run_HSM_analysis.bat` and `run_HSM_analysis.ps1` also contain the original
   absolute local workspace path and must be replaced with portable launchers.
5. Verify that the Python baseline drivers reproduce the reported GTT, BS, IGES,
   and saiyanH result formats after a clean environment setup. Document external
   software that must be installed separately.
6. Export the final all-method Wilcoxon test outputs based on the ten runs used
   in the revised paper and place them beside the raw result CSV files.
7. Confirm the redistribution terms and citations for the eight benchmark
   datasets. If redistribution is not permitted, provide a download and
   preparation script instead of the files themselves.

## Verification status

- The HSM cache-size heatmap script runs successfully from
  `plotting/hsm_cache_size/` using the included CSV data.
- The core plotting entry point requires `matplotlib` and the other packages
  listed in `plotting/requirements.txt`; those packages are not installed in
  the current verification runtime.
- MATLAB R2021a is available on the preparation machine, but BNT is not copied
  into this staging package and has not yet been tested from a clean checkout.

## Included final-result material

- The 24 benchmark `.mat` datasets for eight networks at sample sizes 500, 1000,
  and 3000.
- Original paper comparison, ablation, convergence, transferability, and prior
  statistical-test results.
- Newly added GTT, BS, IGES, and saiyanH result CSV files.
- The latest merged HSM capacity data and the six final cache-size heatmaps.

## Intentionally excluded

- Third-party BNT, pgmpy, Bayesys, and other external package source trees.
- Temporary cache files, IDE settings, Python bytecode, intermediate previews,
  and obsolete HSM plotting variants.
