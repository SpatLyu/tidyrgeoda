# Spatial C(K)luster Analysis by Tree Edge Removal(SKATER)

A wrapper function for
[`rgeoda::skater()`](https://geodacenter.github.io/rgeoda/reference/skater.html).SKATER
forms clusters by spatially partitioning data that has similar values
for features of interest.

## Usage

``` r
st_skater(
  sfj,
  varcol,
  k,
  wt = NULL,
  boundvar = NULL,
  min_bound = 0,
  scale_method = "standardize",
  distance_method = "euclidean",
  seed = 123456789,
  cpu_threads = 6,
  rdist = numeric()
)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- varcol:

  The variable selected to calculate spatial lag, which is a character.

- k:

  The number of clusters.

- wt:

  (optional) The spatial weights object,which can use
  [`st_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_weights.md)
  to construct,default is constructed by `st_weights(sfj,'contiguity')`.

- boundvar:

  (optional) A data frame / tibble with selected bound variable.

- min_bound:

  (optional) A minimum bound value that applies to all clusters.

- scale_method:

  (optional) One of the scaling methods 'raw', 'standardize', 'demean',
  'mad', 'range_standardize', 'range_adjust' to apply on input data.
  Default is 'standardize' (Z-score normalization).

- distance_method:

  (optional) The distance method used to compute the distance between
  observation i and j. Defaults to "euclidean". Options are "euclidean"
  and "manhattan".

- seed:

  (int,optional) The seed for random number generator. Defaults to
  123456789.

- cpu_threads:

  (optional) The number of cpu threads used for parallel computation.

- rdist:

  (optional) The distance matrix (lower triangular matrix, column wise
  storage).

## Value

A names list with names "Clusters", "Total sum of squares",
"Within-cluster sum of squares", "Total within-cluster sum of squares",
and "The ratio of between to total sum of squares".

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
library(sf)
guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
guerry_clusters = st_skater(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),4)
guerry_clusters
#> $Clusters
#>  [1] 3 2 3 1 1 1 2 1 2 1 1 1 2 1 1 3 3 3 2 4 3 1 2 1 2 2 4 1 1 1 1 1 4 3 4 1 2 1
#> [39] 4 3 3 4 2 1 1 1 4 4 2 2 4 2 2 4 2 3 2 2 4 2 3 1 1 1 2 2 1 2 3 4 2 2 2 2 3 2
#> [77] 1 1 1 1 3 3 3 2 2
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 113.20764  43.15129 103.93006  84.62611
#> 
#> $`Total within-cluster sum of squares`
#> [1] 159.0849
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.3156447
#> 
```
