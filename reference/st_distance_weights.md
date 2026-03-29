# Distance-based Spatial Weights

Create a distance-based weights

## Usage

``` r
st_distance_weights(
  sfj,
  unit = "km",
  dist_thres = NULL,
  power = 1,
  is_inverse = FALSE
)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- unit:

  (optional) The unit for calculating spatial distance, can be
  'km'(default) or 'mile'.

- dist_thres:

  (optional) A positive numeric value of distance threshold.

- power:

  (optional) The power (or exponent) of a number indicates how many
  times to use the number in a multiplication.Default is 1.

- is_inverse:

  (optional) FALSE (default) or TRUE, apply inverse on distance value.

## Value

An instance of rgeoda Weight-class.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_distance_weights(guerry)
#> Reference class object of class "Weight"
#> Field "gda_w":
#> An object of class "p_GeoDaWeight"
#> Slot "pointer":
#> <pointer: 0x55e56c2d3aa0>
#> 
#> Field "is_symmetric":
#> [1] TRUE
#> Field "sparsity":
#> [1] 0.04346021
#> Field "min_neighbors":
#> [1] 1
#> Field "max_neighbors":
#> [1] 7
#> Field "num_obs":
#> [1] 85
#> Field "mean_neighbors":
#> [1] 3.694118
#> Field "median_neighbors":
#> [1] 4
#> Field "has_isolates":
#> [1] FALSE
```
