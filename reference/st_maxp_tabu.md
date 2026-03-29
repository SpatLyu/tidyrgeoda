# A tabu-search algorithm to solve the max-p-region problem

A wrapper function for
[`rgeoda::maxp_tabu()`](https://geodacenter.github.io/rgeoda/reference/maxp_tabu.html).The
max-p-region problem is a special case of constrained clustering where a
finite number of geographical areas are aggregated into the maximum
number of regions (max-p-regions), such that each region is
geographically connected and the clusters could maximize internal
homogeneity.

## Usage

``` r
st_maxp_tabu(
  sfj,
  varcol,
  wt = NULL,
  boundvar,
  min_bound,
  tabu_length = 10,
  conv_tabu = 10,
  iterations = 99,
  initial_regions = vector("numeric"),
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

- wt:

  (optional) The spatial weights object,which can use
  [`st_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_weights.md)
  to construct,default is constructed by `st_weights(sfj,'contiguity')`.

- boundvar:

  A numeric vector of selected bounding variable.

- min_bound:

  A minimum value that the sum value of bounding variable int each
  cluster should be greater than.

- tabu_length:

  (optional) The length of a tabu search heuristic of tabu algorithm.
  Defaults to 10.

- conv_tabu:

  (optional): The number of non-improving moves. Defaults to 10.

- iterations:

  (optional) The number of iterations of greedy algorithm. Defaults to
  99.

- initial_regions:

  (optional) The initial regions that the local search starts with.
  Default is empty. means the local search starts with a random process
  to "grow" clusters.

- scale_method:

  (optional) One of the scaling methods 'raw', 'standardize', 'demean',
  'mad', 'range_standardize', 'range_adjust' to apply on input data.
  Default is 'standardize' (Z-score normalization).

- distance_method:

  (optional) The distance method used to compute the distance betwen
  observation i and j. Defaults to "euclidean". Options are "euclidean"
  and "manhattan"

- seed:

  (optional) The seed for random number generator. Defaults to
  123456789.

- cpu_threads:

  (optional) The number of cpu threads used for parallel
  computation.Default is 6.

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
guerry_clusters = st_maxp_tabu(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
'Infants','Suicids'),boundvar = 'Pop1831',min_bound = 3236.67)
guerry_clusters
#> $Clusters
#>  [1] 2 5 1 3 3 1 5 1 4 1 1 3 6 1 2 7 2 7 5 8 2 7 4 3 6 6 8 3 1 1 7 3 8 2 7 3 5 7
#> [39] 2 2 1 8 6 1 7 3 2 2 4 5 8 4 4 8 4 2 5 5 2 6 1 1 1 1 4 4 3 4 2 2 6 6 5 6 2 6
#> [77] 1 1 3 3 8 7 7 4 5
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 66.19164 11.66559 19.22806 19.64391 53.79476 32.85303 37.57316 22.12871
#> 
#> $`Total within-cluster sum of squares`
#> [1] 240.9211
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.4780181
#> 
```
