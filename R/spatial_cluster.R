#' @title Spatial C(K)luster Analysis by Tree Edge Removal(SKATER)
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::skater()`.SKATER forms clusters by spatially
#' partitioning data that has similar values for features of interest.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance
#' between observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan".
#' @param seed (int,optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda skater
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_skater(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),4)
#' guerry_clusters
st_skater = \(sfj,varcol,k,wt = NULL,boundvar = NULL,min_bound = 0,scale_method = "standardize",
              distance_method = "euclidean",seed = 123456789,cpu_threads = 6,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::skater(k,wt,df,boundvar,min_bound,scale_method,
                        distance_method,seed,cpu_threads,rdist))
}

#' @title Regionalization with dynamically constrained agglomerative clustering and partitioning(REDCAP)
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::redcap()`.REDCAP (Regionalization with dynamically constrained agglomerative clustering and partitioning)
#' is developed by D. Guo (2008). Like SKATER, REDCAP starts from building a spanning tree with 4 different ways (single-linkage, average-linkage,
#' ward-linkage and the complete-linkage). The single-linkage way leads to build a minimum spanning tree. Then,REDCAP provides 2 different ways
#' (first-order and full-order constraining) to prune the tree to find clusters. The first-order approach with a minimum spanning tree is exactly
#' the same with SKATER. In GeoDa and pygeoda, the following methods are provided: \* First-order and Single-linkage \* Full-order and Complete-linkage
#'  \* Full-order and Average-linkage \* Full-order and Single-linkage \* Full-order and Ward-linkage.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param method (optional) "firstorder-singlelinkage", "fullorder-completelinkage",
#' "fullorder-averagelinkage"(default),"fullorder-singlelinkage", "fullorder-wardlinkage"
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance
#' between observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan".
#' @param seed (int,optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda redcap
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_redcap(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),
#' 4,method = "fullorder-completelinkage")
#' guerry_clusters
st_redcap = \(sfj,varcol,k,wt = NULL,boundvar = NULL,method = "fullorder-averagelinkage",min_bound = 0,scale_method = "standardize",
              distance_method = "euclidean",seed = 123456789,cpu_threads = 6,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::redcap(k,wt,df,method,boundvar,min_bound,scale_method,
                        distance_method,seed,cpu_threads,rdist))
}

#' @title Spatially Constrained Hierarchical Clucstering (SCHC)
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::schc()`.Spatially constrained hierarchical clustering is
#' a special form of constrained clustering, where the constraint is based on contiguity (common borders).
#' The method builds up the clusters using agglomerative hierarchical clustering methods: single linkage,
#' complete linkage, average linkage and Ward's method (a special form of centroid linkage). Meanwhile,
#' it also maintains the spatial contiguity when merging two clusters.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param method (optional) "single", "complete", "average"(default),"ward".
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance
#' between observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan".
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda schc
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_schc(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants','Suicids'),
#' 4,method = "complete")
#' guerry_clusters
st_schc = \(sfj,varcol,k,wt = NULL,boundvar = NULL,method = "average",min_bound = 0,
            scale_method = "standardize",distance_method = "euclidean",rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::schc(k,wt,df,method,boundvar,min_bound,
                      scale_method,distance_method,rdist))
}

#' @title A greedy algorithm to solve the AZP problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::azp_greedy()`.The automatic zoning procedure (AZP) was
#' initially outlined in Openshaw (1977) as a way to address some of the consequences of
#' the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to
#' find the best set of combinations of contiguous spatial units into p regions, minimizing
#' the within sum of squares as a criterion of homogeneity. The number of regions needs to
#' be specified beforehand.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL
#' "automatic regionalization with initial seed location".
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda azp_greedy
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_azp_greedy(guerry,c('Crm_prs','Crm_prp','Litercy',
#' 'Donatns','Infants','Suicids'),5)
#' guerry_clusters
st_azp_greedy = \(sfj,varcol,k,wt = NULL,boundvar = NULL,min_bound = 0,inits = 0,initial_regions = vector("numeric"),
                  scale_method = "standardize",distance_method = "euclidean",seed = 123456789,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::azp_greedy(k,wt,df,boundvar,min_bound,inits,initial_regions,
                            scale_method,distance_method,seed,rdist))
}

#' @title A simulated annealing algorithm to solve the AZP problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::azp_sa()`.The automatic zoning procedure (AZP) was
#' initially outlined in Openshaw (1977) as a way to address some of the consequences of
#' the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to
#' find the best set of combinations of contiguous spatial units into p regions, minimizing
#' the within sum of squares as a criterion of homogeneity. The number of regions needs to
#' be specified beforehand.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param cooling_rate (optional) The cooling rate of a simulated annealing algorithm. Defaults to 0.85.
#' @param sa_maxit (optional) The number of iterations of simulated annealing. Defaults to 1.
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL
#' "automatic regionalization with initial seed location".
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda azp_sa
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_azp_sa(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
#' 'Infants','Suicids'),5)
#' guerry_clusters
st_azp_sa = \(sfj,varcol,k,wt = NULL,boundvar = NULL,cooling_rate = 0.85,sa_maxit = 1,min_bound = 0,inits = 0,
              initial_regions = vector("numeric"),scale_method = "standardize",distance_method = "euclidean",seed = 123456789,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::azp_sa(k,wt,df,cooling_rate,sa_maxit,boundvar,min_bound,inits,
                        initial_regions,scale_method,distance_method,seed,rdist))
}

#' @title A tabu algorithm to solve the AZP problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::azp_tabu()`.The automatic zoning procedure (AZP) was
#' initially outlined in Openshaw (1977) as a way to address some of the consequences of
#' the modifiable areal unit problem (MAUP). In essence, it consists of a heuristic to
#' find the best set of combinations of contiguous spatial units into p regions, minimizing
#' the within sum of squares as a criterion of homogeneity. The number of regions needs to
#' be specified beforehand.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k The number of clusters.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar (optional) A data frame / tibble with selected bound variable.
#' @param tabu_length (optional) The length of a tabu search heuristic of tabu algorithm. Defaults to 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param min_bound (optional) A minimum bound value that applies to all clusters.
#' @param inits (optional) The number of construction re-runs, which is for ARiSeL
#' "automatic regionalization with initial seed location".
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda azp_tabu
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_azp_tabu(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
#' 'Infants','Suicids'),5)
#' guerry_clusters
st_azp_tabu = \(sfj,varcol,k,wt = NULL,boundvar = NULL,tabu_length = 10,conv_tabu = 10,min_bound = 0,inits = 0,initial_regions = vector("numeric"),
                scale_method = "standardize",distance_method = "euclidean",seed = 123456789,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }
  if (is.null(boundvar)) {
    boundvar = data.frame()
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))

  return(rgeoda::azp_tabu(k,wt,df,tabu_length,conv_tabu,boundvar,min_bound,inits,
                          initial_regions,scale_method,distance_method,seed,rdist))
}

#' @title A greedy algorithm to solve the max-p-region problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::maxp_greedy()`.The max-p-region problem is a special case
#' of constrained clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is geographically
#' connected and the clusters could maximize internal homogeneity.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar A numeric vector of selected bounding variable.
#' @param min_bound A minimum value that the sum value of bounding variable int each
#' cluster should be greater than.
#' @param iterations (optional) The number of iterations of greedy algorithm. Defaults to 99.
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation.Default is 6.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda maxp_greedy
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_maxp_greedy(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
#' 'Infants','Suicids'),boundvar = 'Pop1831',min_bound = 3236.67)
#' guerry_clusters
st_maxp_greedy = \(sfj,varcol,wt = NULL,boundvar,min_bound,iterations = 99,initial_regions = vector("numeric"),
                   scale_method = "standardize",distance_method = "euclidean",seed = 123456789,cpu_threads = 6,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  boundvar = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(boundvar))

  return(rgeoda::maxp_greedy(wt,df,boundvar,min_bound,iterations,initial_regions,
                             scale_method,distance_method,seed,cpu_threads,rdist))
}

#' @title A simulated annealing algorithm to solve the max-p-region problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::maxp_sa()`.The max-p-region problem is a special case
#' of constrained clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is geographically
#' connected and the clusters could maximize internal homogeneity.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar A numeric vector of selected bounding variable.
#' @param min_bound A minimum value that the sum value of bounding variable int each
#' cluster should be greater than.
#' @param cooling_rate (optional) The cooling rate of a simulated annealing algorithm. Defaults to 0.85.
#' @param sa_maxit (optional) The number of iterations of simulated annealing. Defaults to 1.
#' @param iterations (optional) The number of iterations of greedy algorithm. Defaults to 99.
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation.Default is 6.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda maxp_sa
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_maxp_sa(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
#' 'Infants','Suicids'),boundvar = 'Pop1831',min_bound = 3236.67)
#' guerry_clusters
st_maxp_sa = \(sfj,varcol,wt = NULL,boundvar,min_bound,cooling_rate = 0.85,sa_maxit = 1,iterations = 99,initial_regions = vector("numeric"),
               scale_method = "standardize",distance_method = "euclidean",seed = 123456789,cpu_threads = 6,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  boundvar = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(boundvar))

  return(rgeoda::maxp_sa(wt,df,boundvar,min_bound,cooling_rate,sa_maxit,iterations,
                         initial_regions,scale_method,distance_method,seed,cpu_threads,rdist))
}

#' @title A tabu-search algorithm to solve the max-p-region problem
#' @author Wenbo Lv
#' @description
#' A wrapper function for `rgeoda::maxp_tabu()`.The max-p-region problem is a special case
#' of constrained clustering where a finite number of geographical areas are aggregated into
#' the maximum number of regions (max-p-regions), such that each region is geographically
#' connected and the clusters could maximize internal homogeneity.
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param boundvar A numeric vector of selected bounding variable.
#' @param min_bound A minimum value that the sum value of bounding variable int each
#' cluster should be greater than.
#' @param tabu_length (optional) The length of a tabu search heuristic of tabu algorithm. Defaults to 10.
#' @param conv_tabu (optional): The number of non-improving moves. Defaults to 10.
#' @param iterations (optional) The number of iterations of greedy algorithm. Defaults to 99.
#' @param initial_regions (optional) The initial regions that the local search starts with.
#' Default is empty. means the local search starts with a random process to "grow" clusters.
#' @param scale_method (optional) One of the scaling methods 'raw', 'standardize',
#' 'demean', 'mad', 'range_standardize', 'range_adjust' to apply on input data.
#' Default is 'standardize' (Z-score normalization).
#' @param distance_method (optional) The distance method used to compute the distance betwen
#' observation i and j. Defaults to "euclidean". Options are "euclidean" and "manhattan"
#' @param seed (optional) The seed for random number generator. Defaults to 123456789.
#' @param cpu_threads (optional) The number of cpu threads used for parallel computation.Default is 6.
#' @param rdist (optional) The distance matrix (lower triangular matrix, column wise storage).
#'
#' @return A names list with names "Clusters", "Total sum of squares", "Within-cluster sum of squares",
#' "Total within-cluster sum of squares", and "The ratio of between to total sum of squares".
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda maxp_tabu
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' guerry_clusters = st_maxp_tabu(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns',
#' 'Infants','Suicids'),boundvar = 'Pop1831',min_bound = 3236.67)
#' guerry_clusters
st_maxp_tabu = \(sfj,varcol,wt = NULL,boundvar,min_bound,tabu_length = 10,conv_tabu = 10,iterations = 99,initial_regions = vector("numeric"),
                 scale_method = "standardize",distance_method = "euclidean",seed = 123456789,cpu_threads = 6,rdist = numeric()){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  boundvar = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(boundvar))

  return(rgeoda::maxp_tabu(wt,df,boundvar,min_bound,tabu_length,conv_tabu,iterations,
                          initial_regions,scale_method,distance_method,seed,cpu_threads,rdist))
}
