# Univariate Spatial Stratification

Univariate Spatial Stratification by invoking rgeoda's \*\_breaks
function.

## Usage

``` r
st_breaks(sfj, varcol, break_method = "stddev", k = 6)
```

## Arguments

- sfj:

  An sf, tibble or data.frame object

- varcol:

  The variables selected to run univariate spatial stratification.

- break_method:

  (optional) Which has to be one of "stddev"(default), "hinge15",
  "hinge30", "percentile", "natural", "quantile". When the
  `break_method` is "natural" or "quantile", you can assign the `k`
  argument to have how many breaks, other methods will have 6 groups.

- k:

  (optional)A numeric value indicates how many breaks,default is 6.

## Value

A vector of numeric values of computed breaks

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
library(sf)
guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
st_breaks(guerry,'Crm_prs',break_method = "quantile", k = 5)
#> [1] 13270.5 17704.5 21476.5 26743.5
```
