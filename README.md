# CS-CGA

This repository provides the implementation of CS-CGA, the standard benchmark
BN definitions, the benchmark datasets used in the experiments, and the
reported result files.

| Directory | Contents |
| --- | --- |
| `method/` | MATLAB implementation of CS-CGA, including PCR, EEAM, HSM, benchmark BN constructors, and data-generation code. |
| `data_and_results/` | The 24 benchmark datasets used in the experiments and the reported result files. Newly generated datasets are written here locally and ignored by Git. |

The package intentionally excludes implementations of comparison algorithms,
ablation and supplementary-analysis scripts, and figure-generation scripts.

Before running the code, install the Bayes Net Toolbox (BNT) under
`third_party/bnt-master`. BNT is available at
https://www.cs.ubc.ca/~murphyk/Software/BNT/bnt.html. Then run
`method/matlab/run_CS_CGA_experiments.m`. By default, the script uses the
included benchmark datasets. Set `FlagNewdata` to `true` to generate ten new
datasets for each benchmark BN and sample size.
