#' Create y_list and coords objects for multi_bipps
#'
#' @param x vector of x coordinates
#' @param y vector of y coordinates
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
    dplyr::ungroup()

  max_dim <- c(max(df$X),max(df$Y))
  pixel_width <- max_dim[1]/nx
  pixel_height <- max_dim[2]/ny
    # {
    #   if(is.null(scaling)) {
    #     dplyr::mutate(.,
    #                   X_max = max(X),
    #                   Y_max = max(Y),
    #                   X = X/max(X_max,Y_max),
    #                   Y = Y/max(X_max,Y_max))
    #   } else {
    #     dplyr::mutate(.,
    #                   X = X/scaling_factor,
    #                   Y = Y/scaling_factor)
    #   }
    # }

  total_grid <- expand.grid(gridded_x=1:nx,
                            gridded_y=1:ny) %>%
                bind_cols(expand.grid(x=seq(pixel_width / 2, max_dim[1] - pixel_width / 2, length.out = nx),
                            y=seq(pixel_height / 2,max_dim[2] - pixel_height / 2, length.out = ny)))

  counts <- df %>%
    dplyr::mutate(gridded_x = floor((nx-1)*X/max_dim[1]) + 1,
           gridded_y = floor((ny-1)*Y/max_dim[2]) + 1) %>%
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
    dplyr::distinct(x,y) %>%
    dplyr::select(x,y)

  in_hull <- lapply(unique(image_ids),\(id) {
    hull <- spatstat.geom::convexhull.xy(df %>% dplyr::filter(image_id == id) %>% dplyr::select(X,Y))
    data.frame(in_hull=spatstat.geom::inside.owin(coords$x,coords$y,hull))
  }) %>%
    dplyr::bind_rows()

  counts <- counts %>%
    dplyr::bind_cols(in_hull) %>%
    dplyr::mutate(dplyr::across(-c(image_id,gridded_x,gridded_y,in_hull,x,y),~ifelse(is.na(.x) & in_hull,0,.x))) %>%
    dplyr::select(-in_hull)

  y_list <- counts %>%
    dplyr::group_by(image_id) %>%
    dplyr::group_map(~{
      .x %>%
        dplyr::select(-c(gridded_x,gridded_y,x,y)) %>%
        as.matrix()
    })

  names(y_list) <- unique(image_ids)

  list(y_list=y_list,coords=coords)
}
