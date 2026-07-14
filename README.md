# CS-CGA

This repository provides the implementation of CS-CGA, the standard benchmark
BN definitions, and the result files used in the paper.

| Directory | Contents |
| --- | --- |
| `method/` | MATLAB implementation of CS-CGA, including PCR, EEAM, HSM, benchmark BN constructors, and data-generation code. |
| `data_and_results/` | Experimental result files reported in the paper. Generated benchmark datasets are written here locally and ignored by Git. |

The package intentionally excludes implementations of comparison algorithms,
ablation and supplementary-analysis scripts, and figure-generation scripts.

Before running the code, install the Bayes Net Toolbox (BNT) under
`third_party/bnt-master`. BNT is available at
https://www.cs.ubc.ca/~murphyk/Software/BNT/bnt.html. Then run
`method/matlab/run_CS_CGA_experiments.m`, which generates ten datasets for
each benchmark BN and sample size before executing CS-CGA.
