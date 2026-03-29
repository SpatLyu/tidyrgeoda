# Local Moran with Empirical Bayes(EB) Rate

Function to apply local Moran with EB Rate statistics. The EB rate is
first computed from "event" and "base" variables, and then used in local
moran statistics.

## Usage

``` r
st_local_moran_eb(
  sfj,
  varcol,
  wt = NULL,
  permutations = 999,
  permutation_method = "complete",
  significance_cutoff = 0.05,
  cpu_threads = 6,
  seed = 123456789
)
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

- permutations:

  (optional) The number of permutations for the LISA computation.

- permutation_method:

  (optional) The permutation method used for the LISA computation.
  Options are 'complete', 'lookup'. Default is 'complete'.

- significance_cutoff:

  (optional) A cutoff value for significance p-values to filter
  not-significant clusters.

- cpu_threads:

  (optional) The number of cpu threads used for parallel LISA
  computation.

- seed:

  (optional) The seed for random number generator.

## Value

A factor vector.

## Author

Wenbo Lv <lyu.geosocial@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
guerry = sf::read_sf(guerry_path)
st_local_moran_eb(guerry,c("hr60", "po60"))
} # }
```
