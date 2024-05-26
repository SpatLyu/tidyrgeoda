#' @title Contiguity Spatial Weights
#' @author Wenbo Lv
#' @description
#' Create a contiguity spatial weights with options of "queen", "order",
#' "include lower order" and "precision threshold"
#'
#' @param sfj An sf (simple feature) object.
#' @param queen (Optional) TRUE (default) or FALSE, TRUE implements Queen
#' Contiguity and FALSE implements Rook Contiguity.
#' @param order (Optional) Order of contiguity, default is 1.
#' @param include_lower_order (Optional)  Whether or not the lower order
#' neighbors should be included in the weights structure,default is False.
#' @param precision_threshold (Optional) The precision of the underlying shape
#' file is insufficient to allow for an exact match of coordinates to determine
#' which polygons are neighbors,default is 0.
#'
#' @return An instance of rgeoda Weight-class.
#' @importFrom rgeoda queen_weights
#' @importFrom rgeoda rook_weights
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' queenw = st_contiguity_weights(guerry,queen = TRUE)
st_contiguity_weights = \(sfj,queen = TRUE,order = 1,
                          include_lower_order = FALSE,
                          precision_threshold = 0){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  cont_w = if (queen) {
    rgeoda::queen_weights(sfj,order,include_lower_order,precision_threshold)
  } else {
    rgeoda::rook_weights(sfj,order,include_lower_order,precision_threshold)
  }

  return(cont_w)
}

#' @title Distance-based Spatial Weights
#' @author Wenbo Lv
#' @description Create a distance-based weights
#'
#' @param sfj An sf (simple feature) object.
#' @param unit (optional) The unit for calculating spatial distance, can be
#' 'km'(default) or 'mile'.
#' @param dist_thres (optional) A positive numeric value of distance threshold.
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.Default is 1.
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value.
#'
#' @return An instance of rgeoda Weight-class.
#' @importFrom sf st_is_longlat
#' @importFrom rgeoda min_distthreshold
#' @importFrom rgeoda distance_weights
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_distance_weights(guerry)
st_distance_weights = \(sfj,unit = 'km',dist_thres = NULL,
                        power = 1,is_inverse = FALSE){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  is_arc = ifelse(sf::st_is_longlat(sfj),TRUE,FALSE)
  is_mile = ifelse(unit=='mile',TRUE,FALSE)
  dist_thres = ifelse(is.null(dist_thres),
                      rgeoda::min_distthreshold(sfj,is_arc,is_mile),
                      dist_thres)
  dist_w = rgeoda::distance_weights(sfj,dist_thres,power,is_inverse,is_arc,is_mile)

  return(dist_w)
}

#' @title K-Nearest Neighbors-based Spatial Weights
#' @author Wenbo Lv
#' @description Create a k-nearest neighbors based spatial weights
#'
#' @param sfj An sf (simple feature) object.
#' @param k A positive integer number for k-nearest neighbors.
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.Default is 1.
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value.
#' @param unit (optional) The unit for calculating spatial distance, can be
#' 'km'(default) or 'mile'.
#'
#' @return An instance of rgeoda Weight-class.
#' @importFrom sf st_is_longlat
#' @importFrom rgeoda knn_weights
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_knn_weights(guerry,3)
st_knn_weights = \(sfj,k,power = 1,is_inverse = FALSE,unit = 'km'){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  is_arc = ifelse(sf::st_is_longlat(sfj),TRUE,FALSE)
  is_mile = ifelse(unit=='mile',TRUE,FALSE)
  knn_w = rgeoda::knn_weights(sfj,k,power,is_inverse,is_arc,is_mile)

  return(knn_w)
}

#' @title Distance-based Kernel Spatial Weights
#' @author Wenbo Lv
#' @description
#' Create a kernel weights by specifying a bandwidth and a kernel method
#'
#' @param sfj An sf (simple feature) object.
#' @param kernel (optional) A string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'(default).
#' @param bandwidth (optional) A positive numeric value of bandwidth.
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.Default is 1.
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#' on the diagonal of weights matrix.
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value.
#' @param unit (optional) The unit for calculating spatial distance, can be
#' 'km'(default) or 'mile'.
#'
#' @return An instance of rgeoda Weight-class.
#' @importFrom sf st_is_longlat
#' @importFrom rgeoda min_distthreshold
#' @importFrom rgeoda kernel_weights
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_kernel_weights(guerry)
st_kernel_weights = \(sfj,kernel = "gaussian",bandwidth = NULL,power = 1,
                      use_kernel_diagonals = FALSE,is_inverse = FALSE,unit = 'km'){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  is_arc = ifelse(sf::st_is_longlat(sfj),TRUE,FALSE)
  is_mile = ifelse(unit=='mile',TRUE,FALSE)
  bandwidth = ifelse(is.null(bandwidth),
                     rgeoda::min_distthreshold(sfj,is_arc,is_mile),
                     bandwidth)
  kernel_w = rgeoda::kernel_weights(sfj,bandwidth,kernel,use_kernel_diagonals,
                                    power,is_inverse,is_arc,is_mile)

  return(kernel_w)
}

#' @title K-NN Kernel Spatial Weights
#' @author Wenbo Lv
#' @description
#' Create a kernel weights by specifying k-nearest neighbors and a kernel method
#'
#' @param sfj An sf (simple feature) object.
#' @param k A positive integer number for k-nearest neighbors.
#' @param kernel (optional) A string value, which has to be one of 'triangular',
#' 'uniform', 'epanechnikov', 'quartic', 'gaussian'(default).
#' @param power (optional) The power (or exponent) of a number indicates how
#' many times to use the number in a multiplication.Default is 1.
#' @param adaptive_bandwidth (optional) TRUE (default) or FALSE: TRUE use
#' adaptive bandwidth calculated using distance of k-nearest neithbors,
#' FALSE use max distance of all observation to their k-nearest neighbors.
#' @param use_kernel_diagonals (optional) FALSE (default) or TRUE, apply kernel
#' on the diagonal of weights matrix.
#' @param is_inverse (optional) FALSE (default) or TRUE, apply inverse on
#' distance value.
#' @param unit (optional) The unit for calculating spatial distance, can be
#' 'km'(default) or 'mile'.
#'
#' @return An instance of rgeoda Weight-class.
#' @importFrom sf st_is_longlat
#' @importFrom rgeoda kernel_knn_weights
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_kernel_knn_weights(guerry,6)
st_kernel_knn_weights = \(sfj,k,kernel = "gaussian",power = 1,adaptive_bandwidth = TRUE,
                          use_kernel_diagonals = FALSE,is_inverse = FALSE,unit = 'km'){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  is_arc = ifelse(sf::st_is_longlat(sfj),TRUE,FALSE)
  is_mile = ifelse(unit=='mile',TRUE,FALSE)
  kknn_w = rgeoda::kernel_knn_weights(sfj,k,kernel,adaptive_bandwidth,
                                      use_kernel_diagonals,power,
                                      is_inverse,is_arc,is_mile)

  return(kknn_w)
}

#' @title Construct Spatial Weights
#' @author Wenbo Lv
#' @description Create a spatial weights
#'
#' @param sfj An sf (simple feature) object.
#' @param weight The method used to create spatial weights,which has to be one
#' of 'contiguity', 'distance', 'knn', 'kernel', 'kernel_knn'.
#' @param ... Other arguments to construct spatial weight, see
#' 'tidyrgeoda::st_contiguity_weights','tidyrgeoda::st_distance_weights',
#' 'tidyrgeoda::st_knn_weights','tidyrgeoda::st_kernel_weights',
#' 'tidyrgeoda::st_kernel_knn_weights'.
#'
#' @return An instance of rgeoda Weight-class.
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_weights(guerry,'kernel_knn',6)
st_weights = \(sfj,weight = NULL,...){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  switch(weight,
         "contiguity" = {
           wt = st_contiguity_weights(sfj,...)
         },
         "distance" = {
           wt = st_distance_weights(sfj,...)
         },
         "knn" = {
           wt = st_knn_weights(sfj,...)
         },
         "kernel" = {
           wt = st_kernel_weights(sfj,...)
         },
         "kernel_knn" = {
           wt = st_kernel_knn_weights(sfj,...)
         }
  )

  return(wt)
}

#' @title Summary of Spatial Weights
#' @description Override the summary() function for spatial weights
#' @author Wenbo Lv
#'
#' @param object A Weight object
#' @param ... summary optional parameters
#' @return A summary description of an instance of Weight-class
#' @importFrom tibble tibble
#' @examples
#' \dontrun{
#' library(sf)
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = read_sf(guerry_path)
#' queen_w = tidyrgeoda::st_weights(guerry,'contiguity')
#' summary(queen_w)
#' }
#' @export
summary.Weight = \(object, ...) {
  gda_w = object
  name = c("number of observations:",
            "is symmetric: ",
            "sparsity:",
            "# min neighbors:",
            "# max neighbors:",
            "# mean neighbors:",
            "# median neighbors:",
            "has isolates:")
  value = c(as.character(gda_w$num_obs),
            as.character(gda_w$is_symmetric),
            as.character(gda_w$sparsity),
            as.character(gda_w$min_neighbors),
            as.character(gda_w$max_neighbors),
            as.character(gda_w$mean_neighbors),
            as.character(gda_w$median_neighbors),
            as.character(gda_w$has_isolates))
  output = tibble::tibble(name, value)

  return(output)
}

#' @title Read A Geoda File(.gal,.gwt,.swm)
#' @author Wenbo Lv
#' @description
#' Create a spatial weights object from a geoda file
#'
#' @param file_path The file paht of the geoda file.
#' @param id_vec (optional),the id_vec is the id values used in the geoda file.
#'
#' @return A weights object
#' @importFrom stringr str_sub
#' @importFrom stringr str_to_lower
#' @importFrom rgeoda read_gal
#' @importFrom rgeoda read_gwt
#' @importFrom rgeoda read_swm
#' @export
read_geoda =\(file_path,id_vec = NULL){
  filetype = stringr::str_sub(file_path,-3,-1) %>%
    stringr::str_to_lower()
  switch (filetype,
    'gal' = {
      if (is.null(id_vec)){
        wt = rgeoda::read_gal(file_path)
      } else {
        wt = rgeoda::read_gal(file_path,id_vec)
      }
    },
    'gwt' = {
      if (is.null(id_vec)){
        wt = rgeoda::read_gwt(file_path)
      } else {
        wt = rgeoda::read_gwt(file_path,id_vec)
      }
    },
    'swm' = {
      if (is.null(id_vec)){
        wt = rgeoda::read_swm(file_path)
      } else {
        wt = rgeoda::read_swm(file_path,id_vec)
      }
    }
  )
  return(wt)
}

#' @title Spatial Lag
#' @author Wenbo Lv
#' @description
#' Compute the spatial lag for idx-th observation using selected variable and
#' spatial weights matrix
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#'
#' @return A numeric vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select
#' @importFrom dplyr all_of
#' @importFrom dplyr pull
#' @importFrom magrittr %>%
#' @importFrom rgeoda spatial_lag
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_lag(guerry,'Pop1831')
st_lag = \(sfj,varcol,wt = NULL){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splag = rgeoda::spatial_lag(wt,df) %>%
    dplyr::pull()

  return(splag)
}
