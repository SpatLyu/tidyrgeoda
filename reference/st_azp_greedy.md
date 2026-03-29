# A greedy algorithm to solve the AZP problem

A wrapper function for
[`rgeoda::azp_greedy()`](https://geodacenter.github.io/rgeoda/reference/azp_greedy.html).The
automatic zoning procedure (AZP) was initially outlined in Openshaw
(1977) as a way to address some of the consequences of the modifiable
areal unit problem (MAUP). In essence, it consists of a heuristic to
find the best set of combinations of contiguous spatial units into p
regions, minimizing the within sum of squares as a criterion of
homogeneity. The number of regions needs to be specified beforehand.

## Usage

``` r
st_azp_greedy(
  sfj,
  varcol,
  k,
  wt = NULL,
  boundvar = NULL,
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
guerry_clusters = st_azp_greedy(guerry,c('Crm_prs','Crm_prp','Litercy',
'Donatns','Infants','Suicids'),5)
guerry_clusters
#> $Clusters
#>  [1] 5 2 3 1 1 1 2 1 2 1 1 1 2 1 4 4 3 4 2 3 3 1 2 1 2 2 3 1 1 1 1 1 3 3 3 1 2 1
#> [39] 2 3 1 3 2 1 1 1 3 3 2 2 3 2 2 3 2 3 2 2 3 2 1 1 1 1 2 2 1 2 3 3 2 2 2 2 4 2
#> [77] 1 1 1 1 4 4 4 2 2
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 124.25983  33.99347   0.00000  76.58057  87.79966
#> 
#> $`Total within-cluster sum of squares`
#> [1] 181.3665
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.3598541
#> 
```
