# K-Nearest Neighbors-based Spatial Weights

Create a k-nearest neighbors based spatial weights

## Usage

``` r
st_knn_weights(sfj, k, power = 1, is_inverse = FALSE, unit = "km")
```

## Arguments

- sfj:

  An sf (simple feature) object.

- k:

  A positive integer number for k-nearest neighbors.

- power:

  (optional) The power (or exponent) of a number indicates how many
  times to use the number in a multiplication.Default is 1.

- is_inverse:

  (optional) FALSE (default) or TRUE, apply inverse on distance value.

- unit:

  (optional) The unit for calculating spatial distance, can be
  'km'(default) or 'mile'.

## Value

An instance of rgeoda Weight-class.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_knn_weights(guerry,3)
#> Reference class object of class "Weight"
#> Field "gda_w":
#> An object of class "p_GeoDaWeight"
#> Slot "pointer":
#> <pointer: 0x55e571a8d880>
#> 
#> Field "is_symmetric":
#> [1] FALSE
#> Field "sparsity":
#> [1] 0.03529412
#> Field "min_neighbors":
#> [1] 3
#> Field "max_neighbors":
#> [1] 3
#> Field "num_obs":
#> [1] 85
#> Field "mean_neighbors":
#> [1] 3
#> Field "median_neighbors":
#> [1] 3
#> Field "has_isolates":
#> [1] FALSE
```
