# Summary of Spatial Weights

Warpping the [`summary()`](https://rdrr.io/r/base/summary.html) function
for spatial weights

## Usage

``` r
st_summary(wt, ...)
```

## Arguments

- wt:

  A Weight object

- ...:

  summary optional parameters

## Value

A summary description of an instance of Weight-class

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
library(sf)
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = read_sf(guerry_path)
queen_w = tidyrgeoda::st_weights(guerry,'contiguity')
st_summary(queen_w)
} # }
```
