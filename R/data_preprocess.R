#' Pixellate a grid with spatial data and count points of each type
#'
#' This function takes spatial coordinates of points, scales them to fit within a
#' grid defined by the dimensions `[0, 1] x [0, 1]`, and pixellates the grid. It
#' counts the number of points of each type in each pixel, considering separate images
#' when image IDs are provided.
#'
#' @param x A numeric vector of x-coordinates of the points.
#' @param y A numeric vector of y-coordinates of the points.
#' @param types A character or factor vector indicating the type of each point.
#' @param image_ids A vector of identifiers indicating which image each point belongs to.
#' @param nx An integer specifying the number of pixels in the x-direction.
#' @param ny An integer specifying the number of pixels in the y-direction.
#'
#' @return A list with two elements:
#'   \item{y_list}{A named list of matrices where each matrix corresponds to a unique image ID,
#'                 with rows representing pixels and columns representing counts of each point type.}
#'   \item{coords}{A data frame containing the coordinates of each pixel center.}
#' @export
pixellate_grid <- function(x, y, types, image_ids, nx, ny) {
  # Determine the max range in either x or y for uniform scaling
  max_range <- max(max(x) - min(x), max(y) - min(y))

  # Scaling x and y to fit in [0,1] x [0,1]
  x_scaled <- (x - min(x)) / max_range
  y_scaled <- (y - min(y)) / max_range

  # Create a tibble for the scaled points
  points_df <- tibble::tibble(x = x_scaled, y = y_scaled, image_id = image_ids)

  # Creating the grid boundaries
  x_bins <- seq(0, max(x_scaled), length.out = nx + 1)
  y_bins <- seq(0, max(y_scaled), length.out = ny + 1)

  # Assigning each point to a pixel
  x_pixel <- cut(x_scaled, breaks = x_bins, labels = FALSE, include.lowest = TRUE)
  y_pixel <- cut(y_scaled, breaks = y_bins, labels = FALSE, include.lowest = TRUE)

  # Create a dataframe with the pixel and image information
  data <- data.frame(image_id = image_ids, x_pixel, y_pixel, types)

  # Counting the number of points of each type in each pixel for each image
  pixel_counts <- data %>%
    dplyr::group_by(image_id, x_pixel, y_pixel, types) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = types, values_from = count, values_fill = list(count = 0)) %>%
    dplyr::arrange(image_id, x_pixel, y_pixel)

  # Filling missing combinations with NA
  full_grid <- expand.grid(
    image_id = unique(image_ids),
    x_pixel = 1:nx,
    y_pixel = 1:ny
  )

  # Merge with the counts, ensuring missing combinations are NA
  final_df <- full_grid %>%
    dplyr::left_join(pixel_counts, by = c("image_id", "x_pixel", "y_pixel")) %>%
    dplyr::arrange(image_id, x_pixel, y_pixel)

  # Adjust: Rows with all zeros should be turned to NA for type counts
  type_columns <- setdiff(names(final_df), c("image_id", "x_pixel", "y_pixel"))
  final_df <- final_df %>%
    dplyr::mutate(across(all_of(type_columns), ~ifelse(rowSums(!is.na(select(final_df, all_of(type_columns)))) == 0, NA, .)))

  # Adding pixel coordinates to the dataframe
  final_df <- final_df %>%
    dplyr::mutate(
      x = (x_bins[x_pixel] + x_bins[x_pixel + 1]) / 2,
      y = (y_bins[y_pixel] + y_bins[y_pixel + 1]) / 2
    ) %>%
    dplyr::select(image_id, x, y, everything(), -x_pixel, -y_pixel)

  # Handling in_hull logic (if needed)
  coords <- final_df %>%
    dplyr::distinct(x, y)

  in_hull <- lapply(unique(image_ids), \(id) {
    hull <- spatstat.geom::convexhull.xy(points_df %>% dplyr::filter(image_id == id) %>% dplyr::select(x, y))
    data.frame(in_hull = spatstat.geom::inside.owin(coords$x, coords$y, hull))
  }) %>%
    dplyr::bind_rows()

  # Final adjustment to counts
  counts <- final_df %>%
    dplyr::bind_cols(in_hull) %>%
    dplyr::mutate(dplyr::across(-c(image_id, x, y, in_hull), ~ifelse(is.na(.x) & in_hull, 0, .x))) %>%
    dplyr::select(-in_hull)

  # Grouping the data by image and converting to a list of matrices
  y_list <- counts %>%
    dplyr::group_by(image_id) %>%
    dplyr::group_map(~{
      .x %>%
        dplyr::select(-c(x, y)) %>%
        as.matrix()
    })

  # Naming the list by unique image IDs
  names(y_list) <- unique(image_ids)

  # Return the result as a list
  list(y_list = y_list, coords = coords)
}

