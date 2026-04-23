# Calculate ICAT (Independent Component Analysis Trace) Index for TLS

Quantifies the spatial spread and linear organisation of cells within a
detected TLS. FastICA is applied to the (x, y) coordinates of TLS cells
to estimate independent components; the mixing matrix is used to
reconstruct the data, and the ICAT index is defined as the normalised
trace-standard-deviation of the reconstructed coordinates.

The index is always non-negative because it measures the average spatial
spread per cell rather than the signed trace of the mixing matrix (which
can be negative due to ICA sign ambiguity). Higher values indicate a
more spatially extended, structured cluster.

## Usage

``` r
calc_icat(patientID, tlsID, ldata = NULL)
```

## Arguments

- patientID:

  Character. Sample name in `ldata`.

- tlsID:

  Numeric or integer. TLS identifier (value of `tls_id_knn` for the TLS
  of interest).

- ldata:

  Named list of data frames, or `NULL` to use the global `ldata` object
  (deprecated; pass explicitly).

## Value

A single non-negative numeric value (the ICAT index), or `NA_real_` if
computation is not possible (fewer than 3 cells, or FastICA did not
converge).

## Details

The ICAT index is computed as follows:

1.  Centre the (x, y) coordinates of TLS cells.

2.  Run `fastICA` with 2 components.

3.  Reconstruct \\\hat{X} = S A^T + \mu\\.

4.  Let \\v_1, v_2\\ be the marginal variances of \\\hat{X}\\.

5.  \$\$\text{ICAT} = 100 \times \frac{\sqrt{v_1 + v_2 + 2\sqrt{v_1
    v_2}}}{\text{nrow}(X)}\$\$

If the requested TLS contains fewer than 3 cells, or FastICA does not
converge, the function returns `NA_real_` with an informative message
rather than throwing an error.

## Examples

``` r
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#> Detected TLS: 2
if (max(ldata[["ToySample"]]$tls_id_knn, na.rm = TRUE) > 0) {
  icat <- calc_icat("ToySample", tlsID = 1, ldata = ldata)
  print(icat)
}
#> [1] 9.511555
```
