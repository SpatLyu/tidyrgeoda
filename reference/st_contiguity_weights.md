# Contiguity Spatial Weights

Create a contiguity spatial weights with options of "queen", "order",
"include lower order" and "precision threshold"

## Usage

``` r
st_contiguity_weights(
  sfj,
  queen = TRUE,
  order = 1,
  include_lower_order = FALSE,
  precision_threshold = 0
)
```

## Arguments

- sfj:

  An sf (simple feature) object.

- queen:

  (Optional) TRUE (default) or FALSE, TRUE implements Queen Contiguity
  and FALSE implements Rook Contiguity.

- order:

  (Optional) Order of contiguity, default is 1.

- include_lower_order:

  (Optional) Whether or not the lower order neighbors should be included
  in the weights structure,default is False.

- precision_threshold:

  (Optional) The precision of the underlying shape file is insufficient
  to allow for an exact match of coordinates to determine which polygons
  are neighbors,default is 0.

## Value

An instance of rgeoda Weight-class.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
queenw = st_contiguity_weights(guerry,queen = TRUE)
```
