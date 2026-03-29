# (Multivariate) Colocation Local Join Count Statistics

Function to apply (multivariate) colocation local Join Count statistics

## Usage

``` r
st_local_multijoincount(
  sfj,
  varcol,
  wt = NULL,
  permutations = 999,
  permutation_method = "complete",
  significance_cutoff = 0.05,
  cpu_threads = 6,
  seed = 123456789
)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- varcol:

  The variable selected to calculate spatial lag, which is a character.

- wt:

  (optional) The spatial weights object,which can use
  [`st_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_weights.md)
  to construct,default is constructed by `st_weights(sfj,'contiguity')`.

- permutations:

  (optional) The number of permutations for the LISA computation.

- permutation_method:

  (optional) The permutation method used for the LISA computation.
  Options are 'complete', 'lookup'. Default is 'complete'.

- significance_cutoff:

  (optional) A cutoff value for significance p-values to filter
  not-significant clusters.

- cpu_threads:

  (optional) The number of cpu threads used for parallel LISA
  computation.

- seed:

  (optional) The seed for random number generator.

## Value

A factor vector.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_local_multijoincount(guerry,c('TopWealth','TopWealth', 'TopLit'))
#>  [1] Not significant Not significant Not significant Not significant
#>  [5] Not significant Not significant Not significant Not significant
#>  [9] Not significant Not significant Not significant Not significant
#> [13] Not significant Not significant Not significant Not significant
#> [17] Not significant Not significant Not significant Not significant
#> [21] Not significant Not significant Not significant Not significant
#> [25] Not significant Not significant Not significant Not significant
#> [29] Not significant Not significant Not significant Not significant
#> [33] Not significant Not significant Not significant Not significant
#> [37] Not significant Not significant Not significant Not significant
#> [41] Not significant Not significant Not significant Not significant
#> [45] Not significant Not significant Not significant Not significant
#> [49] Not significant Not significant Not significant Not significant
#> [53] Not significant Not significant Not significant Not significant
#> [57] Not significant Not significant Not significant Not significant
#> [61] Not significant Not significant Not significant Not significant
#> [65] Not significant Not significant Not significant Not significant
#> [69] Not significant Not significant Not significant Not significant
#> [73] Not significant Not significant Not significant Not significant
#> [77] Not significant Not significant Not significant Not significant
#> [81] Not significant Not significant Not significant Not significant
#> [85] Not significant
#> Levels: Not significant Significant Undefined Isolated
```
