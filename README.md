# CS-CGA

This repository provides the implementation of CS-CGA and the benchmark data
and result files used in the paper.

| Directory | Contents |
| --- | --- |
| `method/` | MATLAB implementation of CS-CGA, including PCR, EEAM, and HSM. |
| `data_and_results/` | The 24 benchmark datasets and the result files reported in the paper. |

The package intentionally excludes implementations of comparison algorithms,
ablation and supplementary-analysis scripts, and figure-generation scripts.

Before running the code, install the Bayes Net Toolbox (BNT) under
`third_party/bnt-master` and consult `MISSING_ITEMS.md` for the remaining
publication checks.
