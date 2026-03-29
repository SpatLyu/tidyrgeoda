# Spatial Lag

Compute the spatial lag for idx-th observation using selected variable
and spatial weights matrix

## Usage

``` r
st_lag(sfj, varcol, wt = NULL)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- varcol:

  The variable selected to calculate spatial lag, which is a character.

- wt:

  (optional) The spatial weights object,which can use
  [`st_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_weights.md)
  to construct,default is constructed by `st_weights(sfj,'contiguity')`.

## Value

A numeric vector.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_lag(guerry,'Pop1831')
#>  [1] 455.2900 480.3267 382.0433 246.3175 335.2400 324.2800 388.2233 285.0133
#>  [9] 327.8340 303.6340 280.4000 304.6633 485.8033 323.8767 358.7420 405.5440
#> [17] 272.0800 358.1317 329.5114 501.6567 325.4467 367.5271 325.7050 283.0200
#> [25] 459.1650 385.4517 516.1950 297.5033 274.1317 326.6983 389.0975 330.1025
#> [33] 485.7033 270.3450 337.8020 323.5117 370.0660 410.4200 306.6400 430.4057
#> [41] 340.8000 444.7000 311.1014 330.7717 359.4967 321.5680 379.3825 459.0550
#> [49] 322.8817 335.1350 501.0900 417.4475 338.0180 535.1025 477.8900 351.9900
#> [57] 570.6400 491.1167 433.1683 766.8200 300.0600 275.5633 389.4733 261.6250
#> [65] 413.7050 469.1000 452.8700 320.3480 348.6914 345.6550 386.0350 455.2267
#> [73] 441.9037 444.1800 377.7480 649.9140 329.1720 343.9400 251.4933 305.0900
#> [81] 420.2650 325.9483 322.2517 380.5617 306.7860
```
