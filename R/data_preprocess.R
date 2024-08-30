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
create_y_list <- function(X,Y,types,image_ids,nx,ny) {
  # bottom-left corner at 0,0, scale by max X and Y
  df <- data.frame(X = x, Y = y, type = types,image_id = image_ids) %>%
    group_by(image_id) %>%
    mutate(X = X - min(X),
           Y = Y - min(Y)) %>%
    ungroup() %>%
    mutate(X = X/max(X),
           Y = Y/max(Y))

  total_grid <- expand_grid(gridded_x=1:nx,gridded_y=1:ny)

  counts <- df %>%
    mutate(gridded_x = floor((nx-1)*X) + 1,
           gridded_y = floor((ny-1)*Y) + 1) %>%
    count(gridded_x,gridded_y,type,image_id) %>%
    pivot_wider(names_from = type, values_from = n) %>%
    replace(is.na(.),0) %>%
    group_by(image_id) %>%
    group_modify(~{
      .x %>%
        right_join(total_grid, by = c("gridded_x","gridded_y"))
    }) %>%
    ungroup() %>%
    arrange(image_id,gridded_x,gridded_y)

  counts %>%
    ggplot(aes(gridded_x,gridded_y,fill=factor(granulocytes,ordered=TRUE))) +
    geom_tile() +
    facet_wrap(~image_id)

  y_list <- counts %>%
    group_by(image_id) %>%
    group_map(~{
      .x %>%
        select(-c(gridded_x,gridded_y)) %>%
        as.matrix()
    })

  names(y_list) <- levels(image_ids)

  coords <- counts %>%
    distinct(gridded_x,gridded_y) %>%
    mutate(x = gridded_x / max(gridded_x),
           y = gridded_y / max(gridded_y)) %>%
    select(x,y)

  list(y_list=y_list,coords=coords)
}
