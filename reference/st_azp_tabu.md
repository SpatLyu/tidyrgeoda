# A tabu algorithm to solve the AZP problem

A wrapper function for
[`rgeoda::azp_tabu()`](https://geodacenter.github.io/rgeoda/reference/azp_tabu.html).The
automatic zoning procedure (AZP) was initially outlined in Openshaw
(1977) as a way to address some of the consequences of the modifiable
areal unit problem (MAUP). In essence, it consists of a heuristic to
find the best set of combinations of contiguous spatial units into p
regions, minimizing the within sum of squares as a criterion of
homogeneity. The number of regions needs to be specified beforehand.

## Usage

``` r
st_azp_tabu(
  sfj,
  varcol,
  k,
  wt = NULL,
  boundvar = NULL,
  tabu_length = 10,
  conv_tabu = 10,
  min_bound = 0,
  inits = 0,
  initial_regions = vector("numeric"),
  scale_method = "standardize",
  distance_method = "euclidean",
  seed = 123456789,
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

- tabu_length:

  (optional) The length of a tabu search heuristic of tabu algorithm.
  Defaults to 10.

- conv_tabu:

  (optional): The number of non-improving moves. Defaults to 10.

- min_bound:

  (optional) A minimum bound value that applies to all clusters.

- inits:

  (optional) The number of construction re-runs, which is for ARiSeL
  "automatic regionalization with initial seed location".

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
guerry_clusters = st_azp_tabu(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
'Infants','Suicids'),5)
guerry_clusters
#> $Clusters
#>  [1] 4 1 2 3 3 3 1 3 1 3 3 3 2 4 2 2 2 2 1 2 4 2 1 3 1 1 2 3 3 3 2 3 2 2 2 3 1 2
#> [39] 2 4 4 2 1 3 2 3 2 2 1 1 2 1 1 2 1 2 1 1 2 1 4 3 3 3 1 1 3 1 4 2 1 1 1 1 5 1
#> [77] 3 3 3 3 5 5 2 1 1
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 79.37910 13.42356 34.96496 99.58938 63.84544
#> 
#> $`Total within-cluster sum of squares`
#> [1] 212.7976
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.4222174
#> 
```
