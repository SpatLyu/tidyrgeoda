# K-NN Kernel Spatial Weights

Create a kernel weights by specifying k-nearest neighbors and a kernel
method

## Usage

``` r
st_kernel_knn_weights(
  sfj,
  k,
  kernel = "gaussian",
  power = 1,
  adaptive_bandwidth = TRUE,
  use_kernel_diagonals = FALSE,
  is_inverse = FALSE,
  unit = "km"
)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- k:

  A positive integer number for k-nearest neighbors.

- kernel:

  (optional) A string value, which has to be one of 'triangular',
  'uniform', 'epanechnikov', 'quartic', 'gaussian'(default).

- power:

  (optional) The power (or exponent) of a number indicates how many
  times to use the number in a multiplication.Default is 1.

- adaptive_bandwidth:

  (optional) TRUE (default) or FALSE: TRUE use adaptive bandwidth
  calculated using distance of k-nearest neithbors, FALSE use max
  distance of all observation to their k-nearest neighbors.

- use_kernel_diagonals:

  (optional) FALSE (default) or TRUE, apply kernel on the diagonal of
  weights matrix.

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
st_kernel_knn_weights(guerry,6)
#> Reference class object of class "Weight"
#> Field "gda_w":
#> An object of class "p_GeoDaWeight"
#> Slot "pointer":
#> <pointer: 0x55e574561ff0>
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
