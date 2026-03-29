# Local Multivariate Geary Statistics

Function to apply local Multivariate Geary statistics

## Usage

``` r
st_local_multigeary(
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

  The variables selected to calculate spatial lag, which is a character
  vector.

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
st_local_multigeary(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants',
'Suicids'))
#>  [1] Not significant Positive        Positive        Positive       
#>  [5] Not significant Positive        Positive        Positive       
#>  [9] Positive        Positive        Positive        Positive       
#> [13] Not significant Not significant Not significant Not significant
#> [17] Not significant Positive        Positive        Positive       
#> [21] Positive        Positive        Not significant Not significant
#> [25] Positive        Positive        Positive        Positive       
#> [29] Not significant Positive        Not significant Not significant
#> [33] Positive        Not significant Not significant Not significant
#> [37] Not significant Not significant Positive        Not significant
#> [41] Positive        Negative        Positive        Positive       
#> [45] Positive        Not significant Not significant Not significant
#> [49] Positive        Positive        Positive        Not significant
#> [53] Positive        Positive        Not significant Not significant
#> [57] Positive        Positive        Not significant Positive       
#> [61] Positive        Not significant Positive        Not significant
#> [65] Positive        Positive        Not significant Positive       
#> [69] Not significant Positive        Positive        Not significant
#> [73] Positive        Positive        Not significant Positive       
#> [77] Positive        Positive        Positive        Positive       
#> [81] Not significant Positive        Not significant Positive       
#> [85] Positive       
#> Levels: Not significant Positive Negative Undefined Isolated
```
