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
