# Spatially Constrained Hierarchical Clucstering (SCHC)

A wrapper function for
[`rgeoda::schc()`](https://geodacenter.github.io/rgeoda/reference/schc.html).Spatially
constrained hierarchical clustering is a special form of constrained
clustering, where the constraint is based on contiguity (common
borders). The method builds up the clusters using agglomerative
hierarchical clustering methods: single linkage, complete linkage,
average linkage and Ward's method (a special form of centroid linkage).
Meanwhile, it also maintains the spatial contiguity when merging two
clusters.

## Usage

``` r
st_schc(
  sfj,
  varcol,
  k,
  wt = NULL,
  boundvar = NULL,
  method = "average",
  min_bound = 0,
  scale_method = "standardize",
  distance_method = "euclidean",
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

  (optional) "single", "complete", "average"(default),"ward".

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
guerry_clusters = st_schc(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),
4,method = "complete")
guerry_clusters
#> $Clusters
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 4 1 1 1 1 1 1 1 3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [77] 1 1 1 1 1 1 1 1 1
#> 
#> $`Total sum of squares`
#> [1] 504
#> 
#> $`Within-cluster sum of squares`
#> [1] 424.2265   0.0000   0.0000   0.0000
#> 
#> $`Total within-cluster sum of squares`
#> [1] 79.77355
#> 
#> $`The ratio of between to total sum of squares`
#> [1] 0.1582809
#> 
```
