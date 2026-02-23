## Dependencies

All necessary functions are now included in the `src/` folder. The package provides:

- `funcSimulaDataGen6.m` – generates the nine simulated data scenarios.
- `funcL2WMON.m` – weighted isotonic regression (PAV) for a segment.
- `funcXTREMA.m` – finds local extrema indices for the DP algorithm.
- `prefixpeakfinder2.m` – cumulative forward isotonic errors (used for turning‑point detection).

No external toolboxes are required beyond standard MATLAB. The code has been tested with R2020b and later.