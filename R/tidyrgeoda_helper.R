#' @title LISA Fill Scales
#' @author Wenbo Lv
#' @description
#' Provide `ggplot2` fill scales like geoda software.Now it achieve by using
#' `?ggplot2::scale_fill_manual()`.Another achieve can see
#' https://stackoverflow.com/questions/43440068/ggplot2-fix-colors-to-factor-l.
#'
#' @param name The name of the LISA fill scales legend,default is `LISA`
#' @param ... Adjust other legend details for the LISA fill scales, like `labels`
#' to adjust the appearance of legend, see `?ggplot2::scale_fill_manual()`.
#'
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' library(sf)
#' library(ggplot2)
#' guerry_path = system.file("extdata", "Guerry.shp", package = "rgeoda")
#' guerry = read_sf(guerry_path)
#' guerry |>
#'  dplyr::mutate(lisa = st_local_moran(guerry,'Crm_prs')) |>
#'  dplyr::select(lisa) |>
#'  ggplot() +
#'  geom_sf(aes(fill = lisa),lwd = .1,color = 'grey') +
#'  scale_fill_lisa()
scale_fill_lisa = \(name = 'LISA',...){
  ggplot2::scale_fill_manual(
    name,...,
    values = stats::setNames(c("#eeeeee","#FF0000","#0000FF","#a7adf9",
                               "#f4ada8","#464646","#999999"),
                             c("Not significant","High-High","Low-Low","Low-High",
                               "High-Low","Undefined","Isolated")))
}

#' @title Univariate Spatial Stratification
#' @author Wenbo Lv
#' @description
#' Univariate Spatial Stratification by invoking rgeoda's *_breaks function.
#'
#' @param sfj An sf, tibble or data.frame object
#' @param varcol The variables selected to run univariate spatial stratification.
#' @param break_method (optional) Which has to be one of "stddev"(default), "hinge15",
#' "hinge30", "percentile", "natural", "quantile". When the `break_method` is "natural"
#' or "quantile", you can assign the `k` argument to have how many breaks, other methods
#' will have 6 groups.
#' @param k (optional)A numeric value indicates how many breaks,default is 6.
#'
#' @return A vector of numeric values of computed breaks
#' @export
#'
#' @examples
#' library(sf)
#' guerry = read_sf(system.file("extdata", "Guerry.shp", package = "rgeoda"))
#' st_breaks(guerry,'Crm_prs',break_method = "quantile", k = 5)
st_breaks = \(sfj,varcol,break_method = 'stddev',k = 6){
  if (inherits(sfj,"sf")) {
    df = sfj %>%
      sf::st_drop_geometry() %>%
      dplyr::select(dplyr::all_of(varcol))
  } else {
    df = sfj %>%
      dplyr::select(dplyr::all_of(varcol))
  }

  switch(break_method,
         "stddev" = {
           bv = rgeoda::stddev_breaks(df)
         },
         "hinge15" = {
           bv = rgeoda::hinge15_breaks(df)
         },
         "hinge30" = {
           bv = rgeoda::hinge30_breaks(df)
         },
         "percentile" = {
           bv = rgeoda::percentile_breaks(df)
         },
         "natural" = {
           bv = rgeoda::natural_breaks(k,df)
         },
         "quantile" = {
           bv = rgeoda::quantile_breaks(k,df)
         }
  )

  return(bv)
}
