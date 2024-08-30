#' Plot each image in y_list, faceted by types
#'
#' @param y_list y_list from create_y_list
#' @param coords coords from create_y_list
#'
#' @return list of ggplot2 plots
#' @export
plot_y_list <- function(y_list,coords) {
  plots <- lapply(y_list,\(y) {
    dplyr::bind_cols(y,coords) %>%
      tidyr::pivot_longer(-c(x,y),names_to = "type",values_to = "count") %>%
      ggplot2::ggplot(ggplot2::aes(x,y,fill=count)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~type) +
      ggplot2::scale_fill_viridis_c(name="magma")
  })
}
