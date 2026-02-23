import os

files = {
    "README.md": """# ACMR – Adjacency‑Constrained Monotonic Regression

This repository contains MATLAB code for the paper *"Adjacency‑constrained monotonic regression for nonparametric time‑series decomposition"* by Amirhossein Sadoghi.

## Contents

- `run_simulation.m` – main script to run all simulations.
- `run_residual_tests.m` – computes nonparametric tests on residuals.
- `src/` – all algorithm and helper functions.

## Dependencies

All necessary functions are now included in the `src/` folder. The package provides:

- `funcSimulaDataGen6.m` – generates the nine simulated data scenarios.
- `funcL2WMON.m` – weighted isotonic regression (PAV) for a segment.
- `funcXTREMA.m` – finds local extrema indices for the DP algorithm.
- `prefixpeakfinder2.m` – cumulative forward isotonic errors (used for turning‑point detection).

No external toolboxes are required beyond standard MATLAB. The code has been tested with R2020b and later.

## Usage

1. Place all files in a MATLAB working directory.
2. Run `run_simulation(N, K, M, pow)` to generate results for a given sample size `N`, number of periods `K`, repetitions `M`, and trend power `pow` (1 for linear, ~0.7 for nonlinear).
3. The output matrix `results` contains the same information as Table 2 in the paper.
4. Use `run_residual_tests` on the outputs to obtain p‑values for the residual tests (Table 1 in the Internet Appendix).

Example:
```matlab
[res, Yraw] = run_simulation(1600, 266, 5, 1);   % linear trend, n = 1600
p = run_residual_tests(res, Yraw, ...);          % requires all fitted values