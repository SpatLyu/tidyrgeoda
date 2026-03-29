# Regionalization with dynamically constrained agglomerative clustering and partitioning(REDCAP)

A wrapper function for
[`rgeoda::redcap()`](https://geodacenter.github.io/rgeoda/reference/redcap.html).REDCAP
(Regionalization with dynamically constrained agglomerative clustering
and partitioning) is developed by D. Guo (2008). Like SKATER, REDCAP
starts from building a spanning tree with 4 different ways
(single-linkage, average-linkage, ward-linkage and the
complete-linkage). The single-linkage way leads to build a minimum
spanning tree. Then,REDCAP provides 2 different ways (first-order and
full-order constraining) to prune the tree to find clusters. The
first-order approach with a minimum spanning tree is exactly the same
with SKATER. In GeoDa and pygeoda, the following methods are provided:
\\ First-order and Single-linkage \\ Full-order and Complete-linkage \\
Full-order and Average-linkage \\ Full-order and Single-linkage \\
Full-order and Ward-linkage.

## Usage

``` r
st_redcap(
  sfj,
  varcol,
  k,
  wt = NULL,
  boundvar = NULL,
  method = "fullorder-averagelinkage",
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

- method:

  (optional) "firstorder-singlelinkage", "fullorder-completelinkage",
  "fullorder-averagelinkage"(default),"fullorder-singlelinkage",
  "fullorder-wardlinkage"

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
guerry_clusters = st_redcap(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),
4,method = "fullorder-completelinkage")
guerry_clusters
#> $Clusters
#>  [1] 1 2 1 3 3 1 2 3 2 3 3 3 2 1 1 1 1 1 2 1 1 1 2 3 2 2 1 3 3 3 3 3 1 1 1 3 2 3
#> [39] 1 1 4 1 2 3 3 3 1 1 2 2 1 2 2 1 2 1 2 2 1 2 1 3 3 3 2 2 3 2 1 1 2 2 2 2 1 2
#> [77] 3 3 3 3 1 1 1 2 2
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 161.61898  84.62611   0.00000  80.80101
#> 
#> $`Total within-cluster sum of squares`
#> [1] 176.9539
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.351099
#> 
```
