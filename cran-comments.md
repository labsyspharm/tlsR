# CRAN submission comments — tlsR 0.2.0

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.2
* win-builder (R-devel)
* R-hub: windows-x86_64-devel, ubuntu-gcc-devel, macos-arm64

## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes since last version (0.1.2)

This is a feature and bug-fix release. Key changes:

* Two new exported functions: `summarize_TLS()` and `plot_TLS()`.
* Phenotype matching now accepts both singular ("B cell") and plural
  ("B cells") label conventions.
* `calc_icat()` returns `NA` gracefully instead of throwing an error when
  a TLS has too few cells.
* `detect_tic()` now returns visibly and exposes HDBSCAN tuning parameters.
* `scan_clustering()` return value extended with metadata; robust against
  small-tissue edge cases.
* `plot_TLS()` axis labels use plain ASCII to avoid locale encoding issues.
* `VignetteBuilder: knitr` added to DESCRIPTION.
* `LICENSE.md` and `cran-comments.md` excluded via `.Rbuildignore`.
* Comprehensive testthat tests and a workflow vignette added.

Full details in NEWS.md.
