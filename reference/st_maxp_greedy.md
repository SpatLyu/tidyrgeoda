# A greedy algorithm to solve the max-p-region problem

A wrapper function for
[`rgeoda::maxp_greedy()`](https://geodacenter.github.io/rgeoda/reference/maxp_greedy.html).The
max-p-region problem is a special case of constrained clustering where a
finite number of geographical areas are aggregated into the maximum
number of regions (max-p-regions), such that each region is
geographically connected and the clusters could maximize internal
homogeneity.

## Usage

``` r
st_maxp_greedy(
  sfj,
  varcol,
  wt = NULL,
  boundvar,
  min_bound,
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
guerry_clusters = st_maxp_greedy(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
'Infants','Suicids'),boundvar = 'Pop1831',min_bound = 3236.67)
guerry_clusters
#> $Clusters
#>  [1] 7 2 7 1 1 1 2 3 6 3 1 1 8 7 4 4 5 7 2 4 7 3 6 1 8 8 4 1 3 3 3 1 4 5 5 7 2 3
#> [39] 5 7 7 4 8 1 3 1 5 5 6 6 5 6 2 4 6 5 2 2 5 8 7 3 3 3 6 6 2 6 2 5 8 8 2 8 4 8
#> [77] 1 1 1 1 4 4 4 6 2
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 15.92174 36.87513 19.90470 19.85472 35.54800 48.16380 50.69142 32.85303
#> 
#> $`Total within-cluster sum of squares`
#> [1] 244.1875
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.4844989
#> 
```
