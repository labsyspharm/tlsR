# tlsR: Detection and Spatial Analysis of Tertiary Lymphoid Structures

Fast, reproducible detection and quantitative analysis of tertiary
lymphoid structures (TLS) in multiplexed tissue imaging data.

## Typical workflow

1.  Load or prepare a named list of data frames (`ldata`), one per
    tissue sample. Each data frame must contain columns `x`, `y`
    (spatial coordinates in microns), and `phenotype` (character:
    `"B cell"` / `"T cell"` / other).

2.  Run
    [`detect_TLS`](https://amiryousefilab.github.io/tlsR/reference/detect_TLS.md)
    to label B+T co-localised regions.

3.  (Optional) Run
    [`scan_clustering`](https://amiryousefilab.github.io/tlsR/reference/scan_clustering.md)
    to identify windows of significant immune clustering via local
    Ripley's L.

4.  Run
    [`calc_icat`](https://amiryousefilab.github.io/tlsR/reference/calc_icat.md)
    to score the internal linearity/organisation of each detected TLS.

5.  Run
    [`detect_tic`](https://amiryousefilab.github.io/tlsR/reference/detect_tic.md)
    to identify T-cell clusters outside TLS.

6.  Use
    [`summarize_TLS`](https://amiryousefilab.github.io/tlsR/reference/summarize_TLS.md)
    to obtain a tidy summary table.

7.  Use
    [`plot_TLS`](https://amiryousefilab.github.io/tlsR/reference/plot_TLS.md)
    to produce publication-ready spatial plots.

## References

Amiryousefi et al. (2025)
[doi:10.1101/2025.09.21.677465](https://doi.org/10.1101/2025.09.21.677465)

## See also

Useful links:

- <https://github.com/labsyspharm/tlsR>

## Author

**Maintainer**: Ali Amiryousefi <ali_amiryousefi@hms.harvard.edu>
([ORCID](https://orcid.org/0000-0002-6317-3860))

Authors:

- Jeremiah Wala <jeremiah_wala@dfci.harvard.edu>
  ([ORCID](https://orcid.org/0000-0001-6591-1620))

Other contributors:

- Peter Sorger ([ORCID](https://orcid.org/0000-0002-3364-1838))
  \[contributor\]
