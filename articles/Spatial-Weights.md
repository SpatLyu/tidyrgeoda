# Spatial Weights in tidyrgeoda

`tidyrgeoda` provides following methods for spatial weights by invoking
`rgeoda`:

- Queen
- Rook
- Distance based
- K-Nearest Neighbor
- Kernel
- Read GAL/GWT/SWM weights

### Load spatial data and nessary r package.

``` r
library(sf)
library(dplyr)
library(tidyrgeoda)

guerry = read_sf(system.file("extdata","Guerry.shp",package = "rgeoda"))
head(guerry)
```

    ## Simple feature collection with 6 features and 29 fields
    ## Geometry type: MULTIPOLYGON
    ## Dimension:     XY
    ## Bounding box:  xmin: 595532 ymin: 1858801 xmax: 975716 ymax: 2564568
    ## Projected CRS: NTF (Paris) / Lambert zone II
    ## # A tibble: 6 × 30
    ##   CODE_DE COUNT AVE_ID_  dept Region Dprtmnt     Crm_prs Crm_prp Litercy Donatns
    ##   <chr>   <dbl>   <dbl> <dbl> <chr>  <chr>         <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 01          1      49     1 E      Ain           28870   15890      37    5098
    ## 2 02          1     812     2 N      Aisne         26226    5521      51    8901
    ## 3 03          1    1418     3 C      Allier        26747    7925      13   10973
    ## 4 04          1    1603     4 E      Basses-Alp…   12935    7289      46    2733
    ## 5 05          1    1802     5 E      Hautes-Alp…   17488    8174      69    6962
    ## 6 07          1    2249     7 S      Ardeche        9474   10263      27    3188
    ## # ℹ 20 more variables: Infants <dbl>, Suicids <dbl>, MainCty <dbl>,
    ## #   Wealth <dbl>, Commerc <dbl>, Clergy <dbl>, Crm_prn <dbl>, Infntcd <dbl>,
    ## #   Dntn_cl <dbl>, Lottery <dbl>, Desertn <dbl>, Instrct <dbl>, Prsttts <dbl>,
    ## #   Distanc <dbl>, Area <dbl>, Pop1831 <dbl>, TopCrm <dbl>, TopLit <dbl>,
    ## #   TopWealth <dbl>, geometry <MULTIPOLYGON [m]>

## Contiguity Weights

Contiguity means that two spatial units share a common border of
non-zero length. Operationally, we can further distinguish between a
rook and a queen criterion of contiguity, in analogy to the moves
allowed for the such-named pieces on a chess board. The queen criterion
is somewhat more encompassing and defines neighbors as spatial units
sharing a common edge or a common vertex.The rook criterion defines
neighbors by the existence of a common edge between two spatial units.

To create a Queen contiguity weights, one can call the function

``` r
st_contiguity_weights(sfj,queen = TRUE)
```

``` r
qw = st_contiguity_weights(guerry)
st_summary(qw)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          TRUE              
    ## 3 "sparsity:"               0.0581314878892734
    ## 4 "min neighbors:"          2                 
    ## 5 "max neighbors:"          8                 
    ## 6 "mean neighbors:"         4.94117647058824  
    ## 7 "median neighbors:"       5                 
    ## 8 "has isolates:"           FALSE

If `queen` is assigned to `FALSE`,`tidyrgeoda` will invoke `rgeoda` to
create a Rook contiguity weights.

``` r
rw = st_contiguity_weights(guerry,queen = F)
st_summary(rw)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          TRUE              
    ## 3 "sparsity:"               0.0581314878892734
    ## 4 "min neighbors:"          2                 
    ## 5 "max neighbors:"          8                 
    ## 6 "mean neighbors:"         4.94117647058824  
    ## 7 "median neighbors:"       5                 
    ## 8 "has isolates:"           FALSE

You can also create high order of contiguity weights by changing `order`
to more than 1.

``` r
rw = st_contiguity_weights(guerry,queen = F,order = 2)
st_summary(rw)
```

    ## # A tibble: 8 × 2
    ##   name                      value            
    ##   <chr>                     <chr>            
    ## 1 "number of observations:" 85               
    ## 2 "is symmetric: "          TRUE             
    ## 3 "sparsity:"               0.104636678200692
    ## 4 "min neighbors:"          2                
    ## 5 "max neighbors:"          14               
    ## 6 "mean neighbors:"         8.89411764705882 
    ## 7 "median neighbors:"       9                
    ## 8 "has isolates:"           FALSE

Create contiguity weights with 1-order and 2-order together:

``` r
rw2 = st_contiguity_weights(guerry,queen = F,order = 2,
                            include_lower_order = T)
st_summary(rw2)
```

    ## # A tibble: 8 × 2
    ##   name                      value            
    ##   <chr>                     <chr>            
    ## 1 "number of observations:" 85               
    ## 2 "is symmetric: "          TRUE             
    ## 3 "sparsity:"               0.162768166089965
    ## 4 "min neighbors:"          4                
    ## 5 "max neighbors:"          21               
    ## 6 "mean neighbors:"         13.8352941176471 
    ## 7 "median neighbors:"       14               
    ## 8 "has isolates:"           FALSE

## Distance Based Weights

The most straightforward spatial weights matrix constructed from a
distance measure is obtained when i and j are considered neighbors
whenever j falls within a critical distance band from i.You can use
[`st_distance_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_distance_weights.md)
to achieve the Distance Based Weights.

``` r
dw = st_distance_weights(guerry)
st_summary(dw)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          TRUE              
    ## 3 "sparsity:"               0.0434602076124567
    ## 4 "min neighbors:"          1                 
    ## 5 "max neighbors:"          7                 
    ## 6 "mean neighbors:"         3.69411764705882  
    ## 7 "median neighbors:"       4                 
    ## 8 "has isolates:"           FALSE

The distance threshold default is generate use
[`rgeoda::min_distthreshold()`](https://geodacenter.github.io/rgeoda/reference/min_distthreshold.html),but
you can assign `dist_thres` argument by hand.

``` r
dw2 = st_distance_weights(guerry,dist_thres = 1.5e5)
st_summary(dw2)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          TRUE              
    ## 3 "sparsity:"               0.0968858131487889
    ## 4 "min neighbors:"          2                 
    ## 5 "max neighbors:"          13                
    ## 6 "mean neighbors:"         8.23529411764706  
    ## 7 "median neighbors:"       9                 
    ## 8 "has isolates:"           FALSE

## K-Nearest Neighbor Weights

A special case of distance based weights is K-Nearest neighbor weights,
in which every observation will have exactly k neighbors. It can be used
to avoid the problem of isolate in distance-band weights when a smaller
cut-off distance is used. To create a KNN weights, we can call the
function
[`st_knn_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_knn_weights.md):

``` r
knn6_w = st_knn_weights(guerry, 6)
st_summary(knn6_w)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          FALSE             
    ## 3 "sparsity:"               0.0705882352941176
    ## 4 "min neighbors:"          6                 
    ## 5 "max neighbors:"          6                 
    ## 6 "mean neighbors:"         6                 
    ## 7 "median neighbors:"       6                 
    ## 8 "has isolates:"           FALSE

## Kernel Weights

Kernel weights apply kernel function to determine the distance decay in
the derived continuous weights kernel. The kernel weights are defined as
a function \\K\_{z}\\ of the ratio between the distance \\d\_{ij}\\ from
\\i\\ to \\j\\, and the bandwidth \\h_i\\, with \\z =
\frac{d\_{ij}}{h_i}\\.

The kernel functions include

- `triangular`
- `uniform`
- `quadratic`
- `epanechnikov`
- `quartic`
- `gaussian`

Two functions are provided in `tidyrgeoda` to create kernel weights.

#### Use `st_kernel_weights()` for Kernel Weights with afixedbandwidth

``` r
kernel_w = st_kernel_weights(guerry,"uniform")
st_summary(kernel_w)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          FALSE             
    ## 3 "sparsity:"               0.0434602076124567
    ## 4 "min neighbors:"          1                 
    ## 5 "max neighbors:"          7                 
    ## 6 "mean neighbors:"         3.69411764705882  
    ## 7 "median neighbors:"       4                 
    ## 8 "has isolates:"           FALSE

#### Use `st_kernel_knn_weights()` for Kernel Weights with adaptive bandwidth

``` r
adptkernel_w = st_kernel_knn_weights(guerry, 6, "uniform")
st_summary(adptkernel_w)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          FALSE             
    ## 3 "sparsity:"               0.0705882352941176
    ## 4 "min neighbors:"          6                 
    ## 5 "max neighbors:"          6                 
    ## 6 "mean neighbors:"         6                 
    ## 7 "median neighbors:"       6                 
    ## 8 "has isolates:"           FALSE

## Create Spatial Weights object by `st_weights()`

[`st_weights()`](https://spatlyu.github.io/tidyrgeoda/reference/st_weights.md)
is a wrapper function for all above `st_*_weights`,you can use it like:

``` r
qw = st_weights(guerry,weight = 'contiguity')
qw
```

    ## Reference class object of class "Weight"
    ## Field "gda_w":
    ## An object of class "p_GeoDaWeight"
    ## Slot "pointer":
    ## <pointer: 0x55a5b9732410>
    ## 
    ## Field "is_symmetric":
    ## [1] TRUE
    ## Field "sparsity":
    ## [1] 0.05813149
    ## Field "min_neighbors":
    ## [1] 2
    ## Field "max_neighbors":
    ## [1] 8
    ## Field "num_obs":
    ## [1] 85
    ## Field "mean_neighbors":
    ## [1] 4.941176
    ## Field "median_neighbors":
    ## [1] 5
    ## Field "has_isolates":
    ## [1] FALSE

``` r
st_summary(qw)
```

    ## # A tibble: 8 × 2
    ##   name                      value             
    ##   <chr>                     <chr>             
    ## 1 "number of observations:" 85                
    ## 2 "is symmetric: "          TRUE              
    ## 3 "sparsity:"               0.0581314878892734
    ## 4 "min neighbors:"          2                 
    ## 5 "max neighbors:"          8                 
    ## 6 "mean neighbors:"         4.94117647058824  
    ## 7 "median neighbors:"       5                 
    ## 8 "has isolates:"           FALSE
