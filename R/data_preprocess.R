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

  grid <- counts %>%
    select(image_id,gridded_x,gridded_y) %>%
    mutate(x = gridded_x / max(gridded_x),
           y = gridded_y / max(gridded_y)) %>%
    group_by(image_id) %>%
    mutate(ix = 1:n())

  list(y_list=y_list,grid=grid)
}

get_mcmc_verbosity <- function(verbose,mcmc_burn,mcmc_thin,mcmc_keep) {
  if(verbose == 0){
    mcmc_print_every <- 0
  } else {
    if(verbose <= 20){
      mcmc_tot <- mcmc_burn + mcmc_thin * mcmc_keep
      mcmc_print_every <- 1+round(mcmc_tot / verbose)
    } else {
      if(is.infinite(verbose)){
        mcmc_print_every <- 1
      } else {
        mcmc_print_every <- verbose
      }
    }
  }
  mcmc_print_every
}

get_family_id <- function(family,q) {
  family <- if(length(family) == 1){ rep(family, q) } else { family }
  family_in <- data.frame(family=family)

  available_families <- data.frame(id=c(1,4), family=c("poisson", "negbinomial")) # CHECK THIS! numeric for family id?

  family_id <- family_in %>%
    left_join(available_families, by=c("family")) %>%
    pull(id)
}

get_n_partition <- function(n_partition,nx,ny) {
  if(!is.null(n_partition)){
    n_partition
  } else {
    sqrt(nx*ny)
  }
}

set_default <- function(value,default) {
  if(is.null(value)) return(default)
  else return(value)
}
