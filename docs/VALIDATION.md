# Packaging validation

Prepared on 11 November 2026.

Local environment: Linux, CPython 3.12, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8 and mpmath 1.3.0.

- Eight unit tests pass, covering folded moments, mixture variance,
  covariance, a nonzero reconstructed quintic, coherent interference,
  missing physical inputs, data checksums and the embedded pure-Python engine.
- The compact reference data reproduce every stored scientific value
  and covariance of the checked variance report exactly.
- The manifest-verified ZIP imports successfully into a clean data folder.
- Both original drivers run their plot-only modes without a new integral.
- Report generation succeeds from the frozen reference JSON.
- Source packaging uses the Git index, excluding the environment,
  local histories, raw data, caches, compiled libraries and Git metadata.

The engine smoke test evaluates one noncollinear `2 -> 2` configuration
at low quadrature order; it is not a new physical production run. Full
integrals were not rerun for packaging. Windows execution is configured
in GitHub Actions and remains to be checked by that hosted job after upload.
