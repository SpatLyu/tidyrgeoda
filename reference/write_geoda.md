# Save Spatial Weights

Save spatial weights to a file

## Usage

``` r
write_geoda(wt, dsn, id_vec = NULL, layer = NULL)
```

## Arguments

- wt:

  A Weight object

- dsn:

  The path of an output weights file

- id_vec:

  (optional) Defines the unique value of each observation when saving a
  weights file. Default is `tibble::tibble(id_v = 1:wt$num_obs)`.

- layer:

  (optional) The name of the layer of input dataset,default is `"tbf"`.

## Value

A boolean value indicates if save successfully or failed

## Author

Wenbo Lv <lyu.geosocial@gmail.com>
