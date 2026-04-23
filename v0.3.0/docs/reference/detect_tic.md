# Detect T-cell Immune Clusters (TIC) on the Tissue

Applies HDBSCAN to T cells that lie outside of previously detected TLS
regions to identify spatially compact T-cell clusters (TIC). Phenotype
labels `"T cell"` and `"T cells"` are both accepted.

## Usage

``` r
detect_tic(sample, min_pts = 10L, min_cluster_size = 10L, ldata = NULL)
```

## Arguments

- sample:

  Character. Sample name in `ldata`.

- min_pts:

  Integer. HDBSCAN `minPts` parameter: minimum cluster size (default
  `10`). Smaller values detect more, smaller clusters.

- min_cluster_size:

  Integer. Minimum number of T cells for a HDBSCAN cluster to be
  retained; smaller clusters are merged back into noise (label `0`).
  Default `10`.

- ldata:

  Named list of data frames, or `NULL` to use the global `ldata` object
  (deprecated; pass explicitly).

## Value

The input `ldata` list with the sample data frame augmented by one new
column:

- `tcell_cluster_hdbscan`:

  Integer. `0` = noise / not a T-cell cluster; positive integer = TIC
  cluster ID. Non-T-cell rows receive `NA`.

## Examples

``` r
data(toy_ldata)
ldata <- detect_TLS("ToySample", k = 30, ldata = toy_ldata)
#> Detected TLS: 2
ldata <- detect_tic("ToySample", ldata = ldata)
#> detect_tic: 178 T-cell cluster(s) detected in 'ToySample'.
table(ldata[["ToySample"]]$tcell_cluster_hdbscan, useNA = "ifany")
#> 
#>      0      1      2      3      4      5      6      7      8      9     10 
#>   3966     13     12     41     15     13     37     44     28     13     31 
#>     11     12     13     14     15     16     17     18     19     20     21 
#>     94    139     18     95     14     48     13     12     15     17     30 
#>     22     23     24     25     26     27     28     29     30     31     32 
#>     27     16     26     25     39     13     14     56     18     72     22 
#>     33     34     35     36     37     38     39     40     41     42     43 
#>     15    117     32     20     17     14     13     29     18     16     14 
#>     44     45     46     47     48     49     50     51     52     53     54 
#>     46     13     23     51     48     89     52     12    377     11     59 
#>     55     56     57     58     59     60     61     62     63     64     65 
#>     10     65    408     46     17     15     44     17     10     20     22 
#>     66     67     68     69     70     71     72     73     74     75     76 
#>     10     11    155     16     10     29     25     16     10     49     24 
#>     77     78     79     80     81     82     83     84     85     86     87 
#>     24     38     47     20     42     10     98     16     16    123     45 
#>     88     89     90     91     92     93     94     95     96     97     98 
#>     16     43     12     17     21     37     91     24     10     50     37 
#>     99    100    101    102    103    104    105    106    107    108    109 
#>     19     56     11     13     11     13     10     50     14     16     16 
#>    110    111    112    113    114    115    116    117    118    119    120 
#>     13     16     23     46     11     90     14     78    126     52     61 
#>    121    122    123    124    125    126    127    128    129    130    131 
#>     16     12     20     34     12     21     17     21     11     50     46 
#>    132    133    134    135    136    137    138    139    140    141    142 
#>     11     33     13     60     35     13     10     19     42     14    118 
#>    143    144    145    146    147    148    149    150    151    152    153 
#>     30     11     22     45     55     16     12     21     23     12     21 
#>    154    155    156    157    158    159    160    161    162    163    164 
#>     19     27     13     17     34     18     34     23     79     89     16 
#>    165    166    167    168    169    170    171    172    173    174    175 
#>     34     58     89     55     10     15     11     10     43     10     54 
#>    176    177    178   <NA> 
#>     20     23     66 312441 
```
