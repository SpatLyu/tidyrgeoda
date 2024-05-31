#' @title Local Moran Statistics
#' @description
#' Function to apply local Moran statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_moran lisa_labels lisa_clusters
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_moran(guerry,'Crm_prp')
st_local_moran = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                   significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_moran(wt,df,permutations,permutation_method,
                               significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Moran with Empirical Bayes(EB) Rate
#' @description
#' Function to apply local Moran with EB Rate statistics. The EB rate is first computed
#' from "event" and "base" variables, and then used in local moran statistics.
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_moran_eb lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' \dontrun{
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_moran_eb(guerry,c("hr60", "po60"))
#' }
st_local_moran_eb = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                      significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_moran_eb(wt,df,permutations,permutation_method,
                                  significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Bivariate Local Moran Statistics
#' @description
#' Function to apply bivariate local Moran statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_bimoran lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_bimoran(guerry,c('Crm_prs','Litercy'))
st_local_bimoran = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                     significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_bimoran(wt,df,permutations,permutation_method,
                                 significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Geary Statistics
#' @description
#' Function to apply local Geary statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_geary lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_geary(guerry,'Crm_prp')
st_local_geary = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                   significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_geary(wt,df,permutations,permutation_method,
                               significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Multivariate Geary Statistics
#' @description
#' Function to apply local Multivariate Geary statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variables selected to calculate spatial lag, which is a character vector.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_multigeary lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_multigeary(guerry,c('Crm_prs','Crm_prp','Litercy','Donatns','Infants',
#' 'Suicids'))
st_local_multigeary = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                        significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_multigeary(wt,df,permutations,permutation_method,
                                    significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Getis-Ord's G Statistics
#' @description
#' Function to apply Getis-Ord's local G statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_g lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_g(guerry,'Crm_prp')
st_local_g = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
               significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_g(wt,df,permutations,permutation_method,
                           significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Getis-Ord's G* Statistics
#' @description
#' Function to apply Getis-Ord's local G*statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_gstar lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_gstar(guerry,'Crm_prp')
st_local_gstar = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                   significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_gstar(wt,df,permutations,permutation_method,
                               significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Local Join Count Statistics
#' @description
#' Function to apply local Join Count statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_joincount lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_joincount(guerry,'Crm_prp')
st_local_joincount = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                       significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_joincount(wt,df,permutations,permutation_method,
                                   significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Bivariate Local Join Count Statistics
#' @description
#' Function to apply local Bivariate Join Count statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_bijoincount lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' library(magrittr)
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path) %>% dplyr::mutate(InvCrm = 1 - TopCrm)
#' st_local_bijoincount(guerry,c("TopCrm", "InvCrm"))
st_local_bijoincount = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                         significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_bijoincount(wt,df,permutations,permutation_method,
                                     significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title (Multivariate) Colocation Local Join Count Statistics
#' @description
#' Function to apply (multivariate) colocation local Join Count statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_multijoincount lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_multijoincount(guerry,c('TopWealth','TopWealth', 'TopLit'))
st_local_multijoincount = \(sfj,varcol,wt = NULL,permutations = 999,permutation_method = "complete",
                            significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_multijoincount(wt,df,permutations,permutation_method,
                                        significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Quantile LISA Statistics
#' @description
#' Function to apply quantile LISA statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k A value indicates the number of quantiles. Value range e.g. `[1, 10]`.
#' @param q A value indicates which quantile or interval used in local join count statistics.
#' Value stars from 1.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_quantilelisa lisa_clusters
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_quantilelisa(guerry,"Crm_prs",k = 4, q = 1)
st_local_quantilelisa = \(sfj,varcol,k,q,wt = NULL,permutations = 999,permutation_method = "complete",
                          significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_quantilelisa(wt,df,k,q,permutations,permutation_method,
                                      significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}

#' @title Multivariate Quantile LISA Statistics
#' @description
#' Function to apply multivariate quantile LISA statistics
#' @author Wenbo Lv \email{lyu.geosocial@gmail.com}
#'
#' @param sfj An sf (simple feature) object.
#' @param varcol The variable selected to calculate spatial lag, which is a character.
#' @param k A value indicates the number of quantiles. Value range e.g. `[1, 10]`.
#' @param q A value indicates which quantile or interval used in local join count statistics.
#' Value stars from 1.
#' @param wt (optional) The spatial weights object,which can use `st_weights()` to
#' construct,default is constructed by `st_weights(sfj,'contiguity')`.
#' @param permutations (optional) The number of permutations for the LISA computation.
#' @param permutation_method (optional) The permutation method used for the LISA computation. Options are 'complete', 'lookup'.
#' Default is 'complete'.
#' @param significance_cutoff (optional) A cutoff value for significance p-values to filter
#' not-significant clusters.
#' @param cpu_threads (optional) The number of cpu threads used for parallel LISA computation.
#' @param seed (optional) The seed for random number generator.
#'
#' @return A factor vector.
#' @importFrom sf st_drop_geometry
#' @importFrom dplyr select all_of
#' @importFrom magrittr %>%
#' @importFrom rgeoda local_multiquantilelisa lisa_clusters lisa_labels
#' @export
#'
#' @examples
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = sf::read_sf(guerry_path)
#' st_local_multiquantilelisa(guerry,c("Crm_prp", "Litercy"),c(4,4), c(1,1))
st_local_multiquantilelisa = \(sfj,varcol,k,q,wt = NULL,permutations = 999,permutation_method = "complete",
                               significance_cutoff = 0.05,cpu_threads = 6,seed = 123456789){
  stopifnot("sfj must be `sf` object" = inherits(sfj,"sf"))

  if (is.null(wt)) {
    wt = st_weights(sfj,'contiguity')
  }

  df = sfj %>%
    sf::st_drop_geometry() %>%
    dplyr::select(dplyr::all_of(varcol))
  splisa = rgeoda::local_multiquantilelisa(wt,df,k,q,permutations,permutation_method,
                                           significance_cutoff,cpu_threads,seed)
  splisalabel = rgeoda::lisa_labels(splisa)
  splisac = rgeoda::lisa_clusters(splisa)
  splisac = factor(splisalabel[splisac + 1],levels = splisalabel)

  return(splisac)
}
