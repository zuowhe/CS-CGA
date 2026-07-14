# Publication Readiness Checklist

1. Add a top-level `LICENSE` after confirming the intended license for the
   CS-CGA code.
2. Add exact instructions for installing BNT under `third_party/bnt-master`
   and running `method/matlab/run_CS_CGA_experiments.m`.
3. Record the MATLAB release and BNT version, then verify the main CS-CGA
   experiment from a clean checkout.
4. Confirm the redistribution terms and citations for the eight benchmark
   networks. If redistribution is not permitted, replace the dataset files with
   a download and preparation script.

## Scope

- Included: CS-CGA source code, standard benchmark BN definitions,
  data-generation code, and reported result files.
- Excluded: comparison-algorithm code, ablation and supplementary-analysis
  code, plotting code, third-party package source trees, temporary cache files,
  IDE settings, and Python bytecode.
