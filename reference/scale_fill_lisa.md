# LISA Fill Scales

Provide `ggplot2` fill scales like geoda software.Now it achieve by
using `?ggplot2::scale_fill_manual()`.Another achieve can see
https://stackoverflow.com/questions/43440068/ggplot2-fix-colors-to-factor-l.

## Usage

``` r
scale_fill_lisa(name = "LISA", ...)
```

## Arguments

- name:

  The name of the LISA fill scales legend,default is `LISA`

- ...:

  Adjust other legend details for the LISA fill scales, like `labels` to
  adjust the appearance of legend, see `?ggplot2::scale_fill_manual()`.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
library(sf)
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
library(ggplot2)
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = read_sf(guerry_path)
guerry |>
 dplyr::mutate(lisa = st_local_moran(guerry,'Crm_prs')) |>
 dplyr::select(lisa) |>
 ggplot() +
 geom_sf(aes(fill = lisa),lwd = .1,color = 'grey') +
 scale_fill_lisa()
```
