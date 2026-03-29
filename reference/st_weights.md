# Construct Spatial Weights

Create a spatial weights

## Usage

``` r
st_weights(sfj, weight = NULL, ...)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- weight:

  The method used to create spatial weights,which has to be one of
  'contiguity', 'distance', 'knn', 'kernel', 'kernel_knn'.

- ...:

  Other arguments to construct spatial weight, see
  'tidyrgeoda::st_contiguity_weights','tidyrgeoda::st_distance_weights',
  'tidyrgeoda::st_knn_weights','tidyrgeoda::st_kernel_weights',
  'tidyrgeoda::st_kernel_knn_weights'.

## Value

An instance of rgeoda Weight-class.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_weights(guerry,'kernel_knn',6)
#> Reference class object of class "Weight"
#> Field "gda_w":
#> An object of class "p_GeoDaWeight"
#> Slot "pointer":
#> <pointer: 0x55e5695d5f40>
#> 
#> Field "is_symmetric":
#> [1] FALSE
#> Field "sparsity":
#> [1] 0.07058824
#> Field "min_neighbors":
#> [1] 6
#> Field "max_neighbors":
#> [1] 6
#> Field "num_obs":
#> [1] 85
#> Field "mean_neighbors":
#> [1] 6
#> Field "median_neighbors":
#> [1] 6
#> Field "has_isolates":
#> [1] FALSE
```
