# Toy Multiplexed Imaging Data

A small synthetic dataset mimicking multiplexed tissue imaging data,
used in package examples and tests. The list contains one sample named
`"ToySample"`.

## Usage

``` r
toy_ldata
```

## Format

A named list with one element:

- `ToySample`:

  A `data.frame` with the following columns:

  `x`

  :   Numeric. X coordinate in microns.

  `y`

  :   Numeric. Y coordinate in microns.

  `phenotype`

  :   Character. Cell phenotype label. Values are `"B cell"`,
      `"T cell"`, and `"Other"`.

## Source

Synthetically generated for package examples.

## References

Amiryousefi et al. (2025)
[doi:10.1101/2025.09.21.677465](https://doi.org/10.1101/2025.09.21.677465)

## Examples

``` r
data(toy_ldata)
str(toy_ldata)
#> List of 1
#>  $ ToySample:'data.frame':   322951 obs. of  4 variables:
#>   ..$ x        : int [1:322951] 423 355 731 814 1415 1847 2623 2626 2625 3433 ...
#>   ..$ y        : int [1:322951] 234 460 38 420 24 54 353 353 357 30 ...
#>   ..$ cflag    : int [1:322951] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ phenotype: chr [1:322951] "Others" "Others" "Others" "Others" ...
table(toy_ldata[["ToySample"]]$phenotype)
#> 
#>           B cells Endothelial cells     Myeloid cells            Others 
#>              4446             15843             35189            150350 
#>     Stromal cells           T cells 
#>            106412             10711 
```
