# Bivariate Local Moran Statistics

Function to apply bivariate local Moran statistics

## Usage

``` r
st_local_bimoran(
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
st_local_bimoran(guerry,c('Crm_prs','Litercy'))
#>  [1] Not significant High-High       High-Low        Not significant
#>  [5] Not significant Not significant High-High       Not significant
#>  [9] Low-High        Not significant Not significant Not significant
#> [13] Not significant Low-Low         Not significant Not significant
#> [17] High-Low        Low-Low         High-High       High-Low       
#> [21] High-Low        High-Low        Low-High        Not significant
#> [25] Not significant Not significant High-Low        Not significant
#> [29] Not significant Not significant Not significant Not significant
#> [33] High-Low        High-Low        Low-Low         Not significant
#> [37] High-High       Not significant Not significant High-Low       
#> [41] Low-Low         Low-Low         Not significant Low-Low        
#> [45] Not significant Not significant High-Low        Not significant
#> [49] Low-High        High-High       Not significant High-High      
#> [53] High-High       High-Low        Low-High        Not significant
#> [57] Not significant Not significant Not significant Not significant
#> [61] Low-Low         Not significant Not significant Not significant
#> [65] Low-High        Low-High        Not significant High-High      
#> [69] Not significant Not significant Not significant Not significant
#> [73] High-High       Low-High        Not significant Not significant
#> [77] Not significant Not significant Not significant Not significant
#> [81] Not significant Low-Low         Low-Low         Low-High       
#> [85] Not significant
#> 7 Levels: Not significant High-High Low-Low Low-High High-Low ... Isolated
```
