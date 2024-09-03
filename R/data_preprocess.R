#' Create y_list and coords objects for multi_bipps
#'
#' @param X vector of x coordinates
#' @param Y vector of y coordinates
#' @param types character vector of types
#' @param image_ids character vector of image IDs
#' @param nx integer, number of pixels on x axis
#' @param ny integer, number of pixels on y axis
#'
#' @return list containg y_list and coords
#' @export
create_y_list <- function(x,y,types,image_ids,nx,ny) {
  # bottom-left corner at 0,0, scale by max X and Y
  df <- data.frame(X = x, Y = y, type = types,image_id = image_ids) %>%
    dplyr::group_by(image_id) %>%
    dplyr::mutate(X = X - min(X),
           Y = Y - min(Y)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(X = X/max(X),
           Y = Y/max(Y))

  total_grid <- expand.grid(gridded_x=1:nx,gridded_y=1:ny)

  counts <- df %>%
    dplyr::mutate(gridded_x = floor((nx-1)*X) + 1,
           gridded_y = floor((ny-1)*Y) + 1) %>%
    dplyr::count(gridded_x,gridded_y,type,image_id) %>%
    tidyr::pivot_wider(names_from = type, values_from = n) %>%
    replace(is.na(.),0) %>%
    dplyr::group_by(image_id) %>%
    dplyr::group_modify(~{
      .x %>%
        dplyr::right_join(total_grid, by = c("gridded_x","gridded_y"))
    }) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(image_id,gridded_x,gridded_y)

  coords <- counts %>%
    dplyr::distinct(gridded_x,gridded_y) %>%
    dplyr::mutate(x = gridded_x / max(gridded_x),
                  y = gridded_y / max(gridded_y)) %>%
    dplyr::select(x,y)

  in_hull <- lapply(unique(image_ids),\(id) {
    hull <- spatstat.geom::convexhull.xy(df %>% filter(image_id == id) %>% select(X,Y))
    tibble(in_hull=inside.owin(coords$x,coords$y,hull))
  }) %>%
    bind_rows()

  counts <- counts %>%
    bind_cols(in_hull) %>%
    mutate(across(-c(image_id,gridded_x,gridded_y,in_hull),~ifelse(is.na(.x) & in_hull,0,.x))) %>%
    select(-in_hull)

  y_list <- counts %>%
    dplyr::group_by(image_id) %>%
    dplyr::group_map(~{
      .x %>%
        dplyr::select(-c(gridded_x,gridded_y)) %>%
        as.matrix()
    })

  names(y_list) <- unique(image_ids)

  list(y_list=y_list,coords=coords)
}
