#' Runs main multi_bipps MCMC sampling
#'
#' @param y_list list of vectors of response variables per image
#' @param x_list list of matrices of covariates per image
#' @param coords matrix of coordinates, common to all images
#' @param k number of latent processes
#' @param family family of the response variables
#' @param axis_partition number of blocks per axis
#' @param block_size number of observations per block
#' @param n_samples number of MCMC samples
#' @param n_burnin number of burn in samples
#' @param n_thin thinning rate
#' @param n_threads number of threads to use
#' @param verbose integer of verbosity level, 0 for none
#' @param settings list of settings for the MCMC
#' @param prior list of prior parameters
#' @param starting list of starting values
#' @param debug list of debug options
#' @param indpart whether to use independent partitioning
#' @param just_preprocess whether to just preprocess the data
#'
#' @return list of MCMC samples
#' @export
multi_bipps <- function(
  x,
  y,
  types,
  image_ids,
  covariates=NULL,
  k=NULL,
  nx = 20,
  ny = 20,
  family = "poisson",
  n_partition = NULL,
  # axis_partition = NULL,
  # block_size = 30, # how many times to partition n_x and n_y by? found such that resulting block size is not larger than 50, but also not too small
  n_samples = 1000,
  n_burnin = 100,
  n_thin = 1,
  n_threads = 4,
  verbose = 0,
  adapting=TRUE,
  ps=TRUE, # remove
  saving=TRUE,
  low_mem=FALSE,
  debug = list(
    hmc=0,
    sample_beta=TRUE,
   sample_tausq=TRUE,
   sample_theta=TRUE,
   sample_w=TRUE,
   sample_lambda=TRUE,
   verbose=FALSE,
   debug=FALSE,
   prior = list(beta=NULL,
                tausq=NULL,
                sigmasq = NULL, # remove
                phi=NULL,
                a=NULL,
                nu = NULL,
                toplim = NULL, # check being used
                btmlim = NULL,
                set_unif_bounds=NULL # check
   ),
   starting = list(beta=NULL,
                   tausq=NULL,
                   theta=NULL,
                   lambda=NULL,
                   v=NULL,
                   a=NULL,
                   nu = NULL,
                   mcmcsd=.05,
                   mcmc_startfrom=0)
   ),
  indpart=FALSE,
  just_preprocess=FALSE
){

  # inputs:
  # x,y,types, image_ids,k=sqrt(n_types) or log(n_types) or maybe larger - which is better in terms of prediction?,list of covariates, n_x=20,n_y=20, (or grid coordinates) check for extra intercept that is passed in and filter out.
  # mean of (cell-level) covariates?
  # covariates in empty space? covariates are aggregated into grids, so should be passed in on at least a finer grid.
  # assume empty is zero.
  # starting and prior should be in debug as well.
  # dealing with non-square grids? empty space?

  # init
  if(verbose > 0){
    cat("Bayesian Meshed GP regression model fit via Markov chain Monte Carlo\n")
  }
  model_tag <- "Bayesian Meshed GP regression model\n
    o --> o --> o
    ^     ^     ^
    |     |     |
    o --> o --> o
    ^     ^     ^
    |     |     |
    o --> o --> o\n(Markov chain Monte Carlo)\n"

  set_default <- function(x, default_val){
    return(
      if(is.null(x)){
        default_val
      } else {
        x
      }
    )
  }

  if(!is.factor(image_ids)){
    image_ids <- as.factor(image_ids)
  }
  if(!is.factor(types)){
    types <- as.factor(types)
  }

  unique_images <- levels(image_ids)
  num_images <- length(unique_images)

  # set MCMC sampling parameter defaults
  # mcmc_keep <- n_samples
  # mcmc_burn <- n_burnin
  # mcmc_thin <- n_thin
  #
  # # no option, just keep at 0
  # # or keep in debug
  # which_hmc    <- settings$hmc %>% set_default(0)
  # if(!(which_hmc %in% c(0,1,2,3,4,6,7))){
  #   warning("Invalid sampling algorithm choice. Choose settings$hmc in {0,1,2,3,4,6,7}. Setting hmc to 0.")
  #   which_hmc <- 0
  # }
  #
  # mcmc_adaptive    <- settings$adapting %>% set_default(TRUE)
  # mcmc_verbose     <- debug$verbose %>% set_default(FALSE)
  # mcmc_debug       <- debug$debug %>% set_default(FALSE)
  # saving <- settings$saving %>% set_default(TRUE)
  # low_mem <- settings$low_mem %>% set_default(FALSE)
  #
  # debugdag <- debug$dag %>% set_default(1)

  df <- data.frame(X = x, Y = y, type = types,image_id = image_ids)

  # bottom-left corner at 0,0
  df <- df %>%
    group_by(image_id) %>%
    mutate(X = X - min(X),
           Y = Y - min(Y)) %>%
    ungroup()

  # find largest x and y dimensions
  max_x <- max(df$X)
  max_y <- max(df$Y)

  df <- mutate(df, X = X/max_x, Y = Y/max_y)

  total_grid <- expand_grid(gridded_x=1:nx,gridded_y=1:ny)

  counts <- df %>%
    mutate(gridded_x = floor((nx-1)*X) + 1,
           gridded_y = floor((ny-1)*Y) + 1) %>%
    count(gridded_x,gridded_y,type,image_id)
    # group_by(type) %>%
    # summarise(max_x = max(gridded_x),
    #           max_y = max(gridded_y),
    #           min_x = min(gridded_x),
    #           min_y = min(gridded_y))
    # pivot_wider(names_from = type, values_from = n)

  # max_grids <- df %>%
  #   group_by(image_id) %>%
  #   summarise(max_x = max(gridded_x),
  #             max_y = max(gridded_y))

  # df %>%
  #   filter(type == "B cells") %>%
  #   ggplot(aes(gridded_x,gridded_y,fill=n)) +
  #   geom_tile() +
  #   facet_grid(type~image_id) +
  #   theme_bw()
  #
  # data.frame(X = x, Y = y, type = types,image_id = image_ids) %>%
  #   filter(type == "B cells") %>%
  #   ggplot(aes(X,Y)) +
  #   geom_point() +
  #   facet_grid(type~image_id) +
  #   theme_bw()


  # find max grids in X and Y for each image
  # for grids outside of those max, set all counts to NA
  # for grids inside of those max, set counts to zero if not present

  counts <- counts %>%
    pivot_wider(names_from = type, values_from = n) %>%
    replace(is.na(.),0) %>%
    group_by(image_id) %>%
    group_modify(~{
      .x %>%
        right_join(total_grid, by = c("gridded_x","gridded_y"))
    }) %>%
    ungroup() %>%
    arrange(image_id,gridded_x,gridded_y)

# df %>%
#     pivot_longer(-c(gridded_x,gridded_y,image_id), names_to = "type", values_to = "n") %>%
#   filter(type == "B cells") %>%
#   # mutate(n = ifelse(n == 0,NA,n)) %>%
#     ggplot(aes(gridded_x,gridded_y,fill=n)) +
#       geom_tile() +
#       facet_grid(type~image_id) +
#       theme_bw()

# df %>%
#   pivot_longer(-c(gridded_x,gridded_y,image_id), names_to = "type", values_to = "n") %>%
#   group_by(type) %>%
#   filter(!is.na(n)) %>%
#   summarise(max_x = max(gridded_x))

  # counts <- lapply(unique_images,\(id) {
  #   df_filt <- df %>%
  #     filter(image_id == id)
  #   window <- owin(
  #     xrange = c(min(df_filt$X), max(df_filt$X)),
  #     yrange = c(min(df_filt$Y), max(df_filt$Y))
  #   )
  #   pp <- ppp(df_filt$X, df_filt$Y, marks = df_filt$type, window = window)
  #   qd <- quadrats(pp, nx = nx,ny = ny)
  #   tl <- tiles(qd)
  #   pixel_area <- mean(tile.areas(qd)) # small numerical differences, remove with mean
  #   cp <- lapply(tl,\(o) gridcentres(o,nx=1,ny=1)) %>%
  #     do.call(rbind,.) %>%
  #     as_tibble(rownames = "tile") %>%
  #     unnest(c(x,y)) %>%
  #     mutate(tile = gsub("Tile ","",tile))
  #
  #   qdc <- quadratcount(split(pp),tess=qd)
  #
  #   qdc <- lapply(qdc,\(q) {
  #     rownames(q) <- paste0("row ",1:ny)
  #     colnames(q) <- paste0("col ",rev(1:nx))
  #     q
  #   })
  #
  #   unique_types <- names(qdc)
  #   lapply(1:length(qdc),\(i) {
  #     o <- qdc[[i]]
  #     nm <- unique_types[i]
  #     as.data.frame(o) %>%
  #       rename_with(~nm,"Freq")
  #     }) %>%
  #     Reduce(\(a,b) inner_join(a,b,by=c("x","y")),.) %>%
  #     unite("tile", y, x, sep = ", ") %>%
  #     inner_join(cp, by = "tile") %>%
  #     mutate(image_id = id) %>%
  #     pivot_longer(-c(x,y,image_id,tile), names_to = "type", values_to = "count") %>%
  #     mutate(area = pixel_area)
  # }) %>%
  #   bind_rows() %>%
  #   mutate(image_id = factor(image_id),
  #          type = factor(type),
  #          tile = factor(tile))


  # remove only zero cell types in a given image - similar to presence-only data, rather than presence-absence. may need to look at simulation.

  # counts %>%
  #   group_by(image_id) %>%
  #   group_map(~{
  #     dat <- .x
  #     dat %>%
  #       ggplot(aes(x, y, fill = count)) +
  #       geom_raster() +
  #       facet_wrap(~type) +
  #       theme_minimal() +
  #       ggtitle(.y) +
  #       scico::scale_fill_scico(palette = "batlow")
  #   })

  # Set up identical square grid for each image (?). Need to scale other parameters to be on non-square grids??
  # window <- owin()
  # shared_qd <- quadrats(window,nx = nx, ny = ny)
  # tl <- tiles(shared_qd)
  # shared_coords <- lapply(tl,\(o) gridcentres(o,nx=1,ny=1)) %>%
  #   do.call(rbind,.) %>%
  #   as_tibble(rownames = "tile") %>%
  #   unnest(c(x,y)) %>%
  #   mutate(tile = gsub("Tile ","",tile)) %>%
  #   rename(x_shared = x, y_shared = y)

  # counts <- counts %>%
  #   left_join(shared_coords, by = "tile")

  y_list <- counts %>%
    # pivot_wider(names_from = type, values_from = count) %>%
    group_by(image_id) %>%
    group_map(~{
      .x %>%
        select(-c(gridded_x,gridded_y)) %>%
        # tibble::column_to_rownames("tile") %>%
        as.matrix()
    })

  stopifnot("Images have differing numbers of outcomes!"=length(unique(sapply(y_list, ncol))) == 1)

  q <- ncol(y_list[[1]])
  p <- length(covariates)
  if(p == 0) p <- 1 # CHECK!
  k <- ifelse(is.null(k), q, k)

  # MCMC verbosity
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

  # family id
  family <- if(length(family) == 1){ rep(family, q) } else { family }
  family_in <- data.frame(family=family)

  available_families <- data.frame(id=c(1,4), family=c("poisson", "negbinomial")) # CHECK THIS! numeric for family id?

  family_id <- family_in %>%
    left_join(available_families, by=c("family")) %>%
    pull(id)

  # for spatial data: matern or expon, for spacetime: gneiting 2002
  #n_par_each_process <- ifelse(dd==2, 1, 3)

  axis_partition <- c(nx/n_partition, ny/n_partition)

  # what are we sampling
  sample_w       <- debug$sample_w %>% set_default(TRUE)
  sample_beta    <- debug$sample_beta %>% set_default(TRUE)
  sample_tausq   <- debug$sample_tausq %>% set_default(TRUE)
  sample_theta   <- debug$sample_theta %>% set_default(TRUE)
  sample_lambda  <- debug$sample_lambda %>% set_default(TRUE)


  # data management pt 2 - messy, need to clean up

  # what to do with different NAs in different images?
  all_na_which <- lapply(y_list, \(y) apply(y, 1, \(i) ifelse(sum(is.na(i))==q,NA,1)))
  simdata_list <- lapply(1:num_images,\(i) {
    data.frame(ix=1:nrow(shared_coords)) %>%
      cbind(shared_coords, y_list[[i]], all_na_which[[i]]) %>%
      as.data.frame()
      # arrange(!!!rlang::syms(paste0("Var", 1:dd)))
  })

  coords <- simdata_list[[1]] %>%
    dplyr::select(dplyr::contains("Var")) %>%
    as.matrix()
  sort_ix     <- simdata_list[[1]]$ix

  fixed_thresholds <- 1:dd %>% lapply(function(i) kthresholdscp(coords[,i], axis_partition[i]))

  # guaranteed to produce blocks using Mv
  system.time(fake_coords_blocking <- coords %>%
                as.matrix() %>%
                gen_fake_coords(fixed_thresholds, 1) )

  # Domain partitioning and gibbs groups

  system.time(coords_blocking_list <- lapply(simdata_list,\(simdata) {
    coords_blocking <- coords %>%
                  as.matrix() %>%
                  tessellation_axis_parallel_fix(fixed_thresholds, 1) %>%
                  dplyr::mutate(na_which = simdata$na_which, ix=sort_ix)

    # check if some blocks come up totally empty
    blocks_prop <- coords_blocking[,paste0("L", 1:dd)] %>% unique()
    blocks_fake <- fake_coords_blocking[,paste0("L", 1:dd)] %>% unique()
    if(nrow(blocks_fake) != nrow(blocks_prop)){
      #cat("Adding fake coords to avoid empty blocks ~ don't like? Use lower [axis_partition]\n")
      # with current Mv, some blocks are completely empty
      # this messes with indexing. so we add completely unobserved coords
      suppressMessages(adding_blocks <- blocks_fake %>% dplyr::setdiff(blocks_prop) %>%
                         dplyr::left_join(fake_coords_blocking))
      coords_blocking <- dplyr::bind_rows(coords_blocking, adding_blocks)

    }
    coords_blocking
  }))

  # DAG
  system.time(parents_children_list <- lapply(coords_blocking_list,\(coords_blocking) {
    if(dd < 4){
      parents_children <- mesh_graph_build(coords_blocking %>% dplyr::select(-ix), axis_partition, FALSE, n_threads, debugdag)
    } else {
      stop("Input dimension is too high?!")
    }
  }))



  parents                      <- lapply(parents_children_list,\(pc) pc[["parents"]])
  children                     <- lapply(parents_children_list,\(pc) pc[["children"]])
  block_names                  <- lapply(parents_children_list,\(pc) pc[["names"]])
  block_groups                 <- lapply(parents_children_list,\(pc) pc[["groups"]])

  # check back on this
  # if(indpart){
  #   parents %<>% lapply(function(x) lapply(x,\(x_i) numeric(0)))
  #   children %<>% lapply(function(x) numeric(0))
  #   block_groups %<>% rep(0, length(block_groups))
  # }

  suppressMessages(simdata_in_list <- lapply(1:length(simdata_list),\(i) {
    simdata_in <- coords_blocking_list[[i]] %>%
      dplyr::select(-na_which) %>%
      dplyr::left_join(simdata_list[[i]])
  }))

  #simdata[is.na(simdata$ix), "ix"] <- seq(nr_start+1, nr_full)


  indexing_list <- lapply(simdata_in_list,\(simdata_in) {
    blocking <- simdata_in$block %>%
      factor() %>% as.integer()
    indexing <- (1:nrow(simdata_in)-1) %>%
      split(blocking)
  })

  # priors and starting values

  if(1){
    # prior and starting values for mcmc

    # nu
    if(is.null(prior$nu)){
      matern_nu <- FALSE
      if(is.null(starting$nu)){
        start_nu <- 0.5
        matern_fix_twonu <- 1
      } else {
        start_nu <- starting$nu
        if(start_nu %in% c(0.5, 1.5, 2.5)){
          matern_fix_twonu <- 2 * start_nu
        }
      }
    } else {
      nu_limits <- prior$nu
      if(length(nu_limits)==1){
        matern_fix_twonu <- floor(nu_limits)*2 + 1
        start_nu <- matern_fix_twonu
        matern_nu <- FALSE
        if(verbose > 0){
          strmessage <- paste0("nu set to ", start_nu/2)
          message(strmessage)
        }
      } else {
        if(diff(nu_limits) == 0){
          matern_fix_twonu <- floor(nu_limits[1])*2 + 1
          start_nu <- matern_fix_twonu
          matern_nu <- F
          if(verbose > 0){
            strmessage <- paste0("nu set to ", start_nu/2)
            message(strmessage)
          }
        } else {
          if(is.null(starting$nu)){
            start_nu <- mean(nu_limits)
          } else {
            start_nu <- starting$nu
          }
          matern_nu <- TRUE
          matern_fix_twonu <- 1 # not applicable
        }
      }

    }

    if(is.null(prior$phi)){
      stop("Need to specify the limits on the Uniform prior for phi via prior$phi.")
    }
    phi_limits <- prior$phi
    if(is.null(starting$phi)){
      start_phi <- mean(phi_limits)
    } else {
      start_phi <- starting$phi
    }

    if(dd == 3){
      if(is.null(prior$a)){
        a_limits <- phi_limits
      } else {
        a_limits <- prior$a
      }
      if(is.null(starting$a)){
        start_a <- mean(a_limits)
      } else {
        start_a <- starting$a
      }
    }

    if(is.null(prior$beta)){
      beta_Vi <- diag(ncol(x_list[[1]])) * 1/100
    } else {
      beta_Vi <- prior$beta
    }

    if(is.null(prior$tausq)){
      tausq_ab <- c(2, 1)
    } else {
      tausq_ab <- prior$tausq
      if(length(tausq_ab) == 1){
        tausq_ab <- c(tausq_ab[1], 0)
      }
    }

    if(is.null(prior$sigmasq)){
      sigmasq_ab <- c(2, 1)
    } else {
      sigmasq_ab <- prior$sigmasq
    }

    btmlim <- prior$btmlim %>% set_default(1e-3)
    toplim <- prior$toplim %>% set_default(1e3)

    # starting values
    if(is.null(starting$beta)){
      start_beta   <- matrix(0, nrow=p, ncol=q)
    } else {
      start_beta   <- starting$beta
    }

    if(is.null(prior$set_unif_bounds)){
      if(dd == 2){
        if(matern_nu){
          theta_names <- c("phi", "nu", "sigmasq")
          npar <- length(theta_names)

          start_theta <- matrix(0, ncol=k, nrow=npar)
          set_unif_bounds <- matrix(0, nrow=npar*k, ncol=2)

          # phi
          set_unif_bounds[seq(1, npar*k, npar),] <- matrix(phi_limits,nrow=1) %x% matrix(1, nrow=k)
          start_theta[1,] <- start_phi

          # nu
          set_unif_bounds[seq(2, npar*k, npar),] <- matrix(nu_limits,nrow=1) %x% matrix(1, nrow=k)
          start_theta[2,] <- start_nu

          # sigmasq expansion
          set_unif_bounds[seq(3, npar*k, npar),] <- matrix(c(btmlim, toplim),nrow=1) %x% matrix(1, nrow=k)
          if((q>1) & (!use_ps)){
            # multivariate without expansion: sigmasq=1 fixed.
            start_theta[3,] <- 1
          } else {
            start_theta[3,] <- btmlim + 1
          }

        } else {
          theta_names <- c("phi", "sigmasq")
          npar <- length(theta_names)

          start_theta <- matrix(0, ncol=k, nrow=npar)

          set_unif_bounds <- matrix(0, nrow=npar*k, ncol=2)
          # phi
          set_unif_bounds[seq(1, npar*k, npar),] <- matrix(phi_limits,nrow=1) %x% matrix(1, nrow=k)
          start_theta[1,] <- start_phi

          # sigmasq expansion
          set_unif_bounds[seq(2, npar*k, npar),] <- matrix(c(btmlim, toplim),nrow=1) %x% matrix(1, nrow=k)
          if((q>1) & (!use_ps)){
            # multivariate without expansion: sigmasq=1 fixed.
            start_theta[2,] <- 1
          } else {
            start_theta[2,] <- btmlim + 1
          }
        }
      } else {
        theta_names <- c("a", "phi", "beta", "sigmasq")
        npar <- length(theta_names)

        start_theta <- matrix(0, ncol=k, nrow=npar)
        set_unif_bounds <- matrix(0, nrow=npar*k, ncol=2)

        # a
        set_unif_bounds[seq(1, npar*k, npar),] <- matrix(a_limits,nrow=1) %x% matrix(1, nrow=k)
        start_theta[1,] <- start_a

        # phi
        set_unif_bounds[seq(2, npar*k, npar),] <- matrix(phi_limits,nrow=1) %x% matrix(1, nrow=k)
        start_theta[2,] <- start_phi

        # beta
        set_unif_bounds[seq(3, npar*k, npar),] <- matrix(c(0,1),nrow=1) %x% matrix(1, nrow=k)
        start_theta[3,] <- 0.5

        # sigmasq expansion
        set_unif_bounds[seq(4, npar*k, npar),] <- matrix(c(btmlim, toplim),nrow=1) %x% matrix(1, nrow=k)
        if((q>1) & (!use_ps)){
          # multivariate without expansion: sigmasq=1 fixed.
          start_theta[4,] <- 1
        } else {
          start_theta[4,] <- btmlim + 1
        }

      }

    } else {
      set_unif_bounds <- prior$set_unif_bounds
    }

    # override defaults if starting values are provided
    if(!is.null(starting$theta)){
      start_theta <- starting$theta
    }

    n_par_each_process <- nrow(start_theta)
    if(is.null(starting$mcmcsd)){
      mcmc_mh_sd <- diag(k * n_par_each_process) * 0.05
    } else {
      if(length(starting$mcmcsd) == 1){
        mcmc_mh_sd <- diag(k * n_par_each_process) * starting$mcmcsd
      } else {
        mcmc_mh_sd <- starting$mcmcsd
      }
    }

    if(is.null(starting$tausq)){
      start_tausq  <- family %>% sapply(function(ff) if(ff == "gaussian"){.1} else {1})
    } else {
      start_tausq  <- starting$tausq
    }

    if(is.null(starting$lambda)){
      start_lambda <- matrix(0, nrow=q, ncol=k)
      diag(start_lambda) <- if(use_ps){10} else {1}
    } else {
      start_lambda <- starting$lambda
    }

    if(is.null(starting$lambda_mask)){
      if(k<=q){
        lambda_mask <- matrix(0, nrow=q, ncol=k)
        lambda_mask[lower.tri(lambda_mask)] <- 1
        diag(lambda_mask) <- 1 #***
      } else {
        stop("starting$lambda_mask needs to be specified")
      }
    } else {
      lambda_mask <- starting$lambda_mask
    }

    if(is.null(starting$mcmc_startfrom)){
      mcmc_startfrom <- 0
    } else {
      mcmc_startfrom <- starting$mcmc_startfrom
    }

    if(is.null(starting$v)){
      start_v <- lapply(1:num_images,\(i) matrix(0, nrow = nrow(simdata_in_list[[1]]), ncol = k))
    } else {
      # this is used to restart MCMC
      # assumes the ordering and the sizing is correct,
      # so no change is necessary and will be input directly to mcmc
      start_v <- starting$v
    }
  }

  # finally prepare data
  y_list <- lapply(1:length(y_list), \(i) {
    y <- simdata_in_list[[i]] %>%
      dplyr::select(dplyr::contains("Y_")) %>%
      as.matrix()
    colnames(y) <- orig_y_colnames
    y
  })


  x_list <- lapply(simdata_in_list,\(simdata_in) {
    x <- simdata_in %>%
      dplyr::select(dplyr::contains("X_")) %>%
      as.matrix()
    colnames(x) <- orig_X_colnames
    x[is.na(x)] <- 0 # NAs if added coords due to empty blocks
    x
  })


  coords <- simdata_in_list[[1]] %>%
    dplyr::select(dplyr::contains("Var")) %>%
    as.matrix()


  coords_renamer <- colnames(coords)
  names(coords_renamer) <- orig_coords_colnames

  coordsdata <- simdata_in_list[[1]] %>%
    dplyr::select(all_of(1:dd), ix) %>%
    dplyr::rename(!!!coords_renamer)

  if(verbose > 0){
    cat("Sending to MCMC.\n")
  }

  mcmc_run <- multi_bipps_mcmc

  if(!just_preprocess){
    comp_time <- system.time({
        results <- mcmc_run(
          y_list,
          family_id,
          x_list,
          coords,
          k,
          parents,
          children,
          block_names,
          block_groups,
          indexing_list,

          set_unif_bounds,
          beta_Vi,


          sigmasq_ab,
          tausq_ab,

          matern_fix_twonu,

          start_v, # this needs to be a list of matrices

          start_lambda,
          lambda_mask,

          start_theta,
          start_beta,
          start_tausq,

          mcmc_mh_sd,

          mcmc_keep, mcmc_burn, mcmc_thin,

          mcmc_startfrom,

          n_threads,

          which_hmc,
          mcmc_adaptive, # adapting

          use_ps,

          mcmc_verbose, mcmc_debug, # verbose, debug
          mcmc_print_every, # print all iter
          low_mem,
          # sampling of:
          # beta tausq sigmasq theta w
          sample_beta, sample_tausq,
          sample_lambda,
          sample_theta, sample_w)
      })
  }

  if(saving){

    listN <- function(...){
      anonList <- list(...)
      names(anonList) <- as.character(substitute(list(...)))[-1]
      anonList
    }

    saved <- listN(y, x, coords_blocking, k,

                   family,
      parents, children,
      block_names, block_groups,

      indexing,

      set_unif_bounds,
      beta_Vi,

      tausq_ab,

      start_v,

      start_lambda,
      lambda_mask,

      start_theta,
      start_beta,
      start_tausq,

      mcmc_mh_sd,

      mcmc_keep, mcmc_burn, mcmc_thin,

      mcmc_startfrom,

      n_threads,

      mcmc_adaptive, # adapting

      use_ps,
      matern_fix_twonu,

      mcmc_verbose, mcmc_debug, # verbose, debug
      mcmc_print_every, # print all iter
      # sampling of:
      # beta tausq sigmasq theta w
      sample_beta, sample_tausq,
      sample_lambda,
      sample_theta, sample_w,
      fixed_thresholds)
  } else {
    saved <- "Model data not saved."
  }

  returning <- list(coordsdata = coordsdata,
                    savedata = saved)


  if(!just_preprocess){
    returning <- c(returning, results)
    class(returning) <- "multi_bipps"
    return ( returning )
  } else {
    class(returning) <- "multi_bipps_preprocess_only"
    return(returning)
  }



  return(returning)

}


