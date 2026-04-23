# Summarize Detected TLS Across Samples

Produces a tidy `data.frame` with one row per sample summarising the
number of detected TLS, their sizes, and (optionally) ICAT scores.

## Usage

``` r
summarize_TLS(ldata, calc_icat_scores = FALSE)
```

## Arguments

- ldata:

  Named list of data frames as returned by
  [`detect_TLS`](https://amiryousefilab.github.io/tlsR/reference/detect_TLS.md)
  (and optionally
  [`detect_tic`](https://amiryousefilab.github.io/tlsR/reference/detect_tic.md)).

- calc_icat_scores:

  Logical. Should ICAT scores be computed for each TLS and appended as a
  list-column? Default `FALSE`.

## Value

A `data.frame` with columns:

- `sample`:

  Sample name.

- `n_TLS`:

  Number of TLS detected.

- `total_cells`:

  Total cells in the sample.

- `TLS_cells`:

  Number of cells assigned to any TLS.

- `TLS_fraction`:

  Fraction of all cells that are TLS cells.

- `mean_TLS_size`:

  Mean cells per TLS (`NA` if n_TLS = 0).

- `n_TIC`:

  Number of T-cell clusters detected by
  [`detect_tic`](https://amiryousefilab.github.io/tlsR/reference/detect_tic.md)
  (`NA` if not yet run).

- `icat_scores`:

  List-column of ICAT scores per TLS (only when
  `calc_icat_scores = TRUE`).

## Examples

``` r
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#> Detected TLS: 2
summarize_TLS(ldata)
#>      sample n_TLS total_cells TLS_cells TLS_fraction mean_TLS_size n_TIC
#> 1 ToySample     2      322951      1568  0.004855226           784    NA
```
