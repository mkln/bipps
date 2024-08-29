#' Runs main multi_bipps MCMC sampling
#'
#' @param x x coordinate of each cell
#' @param y y coordinate of each cell
#' @param types vector of cell types
#' @param image_ids vector of image that each cell comes from
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
  covariates=NULL, # wait on these
  k=NULL,
  nx = 20,
  ny = 20,
  # family = "poisson", # fix to poisson for now
  # n_partition = NULL, # fix to something for now - so that # of grids in each (rectangular) block close to 30-50
  n_samples = 1000,
  n_burnin = 100,
  n_thin = 1, # remove - initially for memory purposes
  n_threads = 4,
  verbose = 0,
  settings = list(adapting=TRUE,
                  saving=TRUE, low_mem=FALSE),
  prior = list(beta=NULL, tausq=NULL, sigmasq = NULL,
               phi=NULL, a=NULL, nu = NULL,
               toplim = NULL, btmlim = NULL, set_unif_bounds=NULL),
  starting = list(beta=NULL, tausq=NULL, theta=NULL,
                  lambda=NULL, v=NULL,  a=NULL, nu = NULL,
                  mcmcsd=.05, mcmc_startfrom=0),
  debug = list(sample_beta=TRUE, sample_tausq=TRUE,
               sample_theta=TRUE, sample_w=TRUE, sample_lambda=TRUE,
               verbose=FALSE, debug=FALSE),
  indpart=FALSE,
  just_preprocess=FALSE
  # adapting=TRUE,
  # saving=TRUE, # restart MCMC later - check what it saves. Still would need to save the last v.
  # low_mem=FALSE, # default to true, or just remove? we really only care about beta, lambda, theta. remove v, w, yhat
  # debug = list(
  #  hmc=0,
  #  sample_beta=TRUE,
  #  sample_tausq=TRUE, # remove for now, since also useful for nb
  #  sample_theta=TRUE,
  #  sample_w=TRUE,
  #  sample_lambda=TRUE,
  #  verbose=FALSE,
  #  debug=FALSE,
  #  prior = list(beta=NULL,
  #               tausq=NULL, # remove, nb
  #               sigmasq = NULL, # remove
  #               phi=NULL, # uniform typically for theta (on ine)
  #               a=NULL, # remove, related to time?
  #               nu = NULL, # matern, remove
  #               toplim = NULL, # check being used, maybe remove
  #               btmlim = NULL,
  #               set_unif_bounds=NULL # check
  #  ),
  #  starting = list(beta=NULL,
  #                  tausq=NULL, # remove
  #                  theta=NULL,
  #                  lambda=NULL,
  #                  v=NULL,
  #                  a=NULL, # remove
  #                  nu = NULL, # remove
  #                  mcmcsd=.05, # used to restart with adaptation
  #                  mcmc_startfrom=0)
  #  ),
  # indpart=FALSE, # put in debug
  # just_preprocess=FALSE # remove, do preprocessing only if mcmc = 0
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

  if(!is.factor(image_ids)){
    image_ids <- as.factor(image_ids)
  }
  if(!is.factor(types)){
    types <- as.factor(types)
  }
  unique_images <- levels(image_ids)
  num_images <- length(unique_images)

  # set MCMC sampling parameter defaults
  #
  # # no option, just keep at 0
  # # or keep in debug
  # which_hmc    <- settings$hmc %>% set_default(0)
  # if(!(which_hmc %in% c(0,1,2,3,4,6,7))){
  #   warning("Invalid sampling algorithm choice. Choose settings$hmc in {0,1,2,3,4,6,7}. Setting hmc to 0.")
  #   which_hmc <- 0
  # }
  #

  mcmc_adaptive    <- settings$adapting %>% set_default(TRUE)
  mcmc_verbose     <- debug$verbose %>% set_default(FALSE)
  mcmc_debug       <- debug$debug %>% set_default(FALSE)
  saving <- settings$saving %>% set_default(TRUE)
  low_mem <- settings$low_mem %>% set_default(FALSE)
  #
  debugdag <- debug$dag %>% set_default(1)

  out <- create_y_list(X,Y,types,image_ids,nx,ny)
  y_list <- out$y_list
  grid <- out$grid

  stopifnot("Images have differing numbers of outcomes!"=length(unique(sapply(y_list, ncol))) == 1)

  q <- ncol(y_list[[1]])
  p <- length(covariates)
  if(p == 0) p <- 1 # CHECK!
  k <- ifelse(is.null(k), q, k)

  # MCMC verbosity
  mcmc_keep <- n_samples
  mcmc_burn <- n_burnin
  mcmc_thin <- n_thin
  mcmc_print_every <- get_mcmc_verbosity(verbose,mcmc_burn,mcmc_thin,mcmc_keep)

  # family id
  family_id <- get_family_id(family,q)

  # for spatial data: matern or expon, for spacetime: gneiting 2002
  #n_par_each_process <- ifelse(dd==2, 1, 3)

  n_partition <- get_n_partition(n_partition, nx, ny)

  axis_partition <- ceiling(c(nx/n_partition, ny/n_partition))

  # what are we sampling - need to fix this
  sample_w       <- debug$sample_w %>% set_default(TRUE)
  sample_beta    <- debug$sample_beta %>% set_default(TRUE)
  sample_tausq   <- debug$sample_tausq %>% set_default(TRUE)
  sample_theta   <- debug$sample_theta %>% set_default(TRUE)
  sample_lambda  <- debug$sample_lambda %>% set_default(TRUE)


  # data management pt 2 - messy, need to clean up

  # what to do with different NAs in different images?
  all_na_which <- lapply(y_list, \(y) apply(y, 1, \(i) ifelse(sum(is.na(i))==q,NA,1)))

  # check if prop of NAs not too large in image
  # fill in NAs in middle of image with zeros

  grid_list <- grid %>%
    group_by(image_id) %>%
    group_split()

  coords <- grid_list[[1]] %>%
    select(x,y)

  # remove, since we are generating the coordinate system
  sort_ix <- grid_list[[1]]$ix


  simdata_list <- lapply(1:num_images,\(i) {
    data.frame(ix=sort_ix) %>%
      cbind(coords, y_list[[i]], all_na_which[[i]]) %>%
      as.data.frame()
  })

  fixed_thresholds <- 1:2 %>% lapply(\(i) kthresholdscp(coords %>% pull(i), axis_partition[i]))

  # guaranteed to produce blocks using Mv - remove, since we already have full grid with NAs
  system.time(fake_coords_blocking <- coords %>%
                as.matrix() %>%
                gen_fake_coords(fixed_thresholds, 1) )

  # Domain partitioning and gibbs groups

  # update this, since only NAs change between images
  system.time(coords_blocking_list <- lapply(all_na_which,\(na_which) {
    coords_blocking <- coords %>%
                  as.matrix() %>%
                  tessellation_axis_parallel_fix(fixed_thresholds, 1) %>%
                  dplyr::mutate(na_which = na_which, ix=sort_ix)

    # check if some blocks come up totally empty
    blocks_prop <- coords_blocking[,paste0("L", 1:2)] %>% unique()
    blocks_fake <- fake_coords_blocking[,paste0("L", 1:2)] %>% unique()
    if(nrow(blocks_fake) != nrow(blocks_prop)){
      # with current Mv, some blocks are completely empty
      # this messes with indexing. so we add completely unobserved coords
      suppressMessages(adding_blocks <- blocks_fake %>%
                         dplyr::setdiff(blocks_prop) %>%
                         dplyr::left_join(fake_coords_blocking))
      coords_blocking <- dplyr::bind_rows(coords_blocking, adding_blocks)

    }
    coords_blocking
  }))

  # DAG
  system.time(parents_children_list <- lapply(coords_blocking_list,\(coords_blocking) {
    parents_children <- mesh_graph_build(coords_blocking %>%
                                           dplyr::select(-ix),
                                         axis_partition,
                                         FALSE,
                                         n_threads,
                                         FALSE)
  }))

  parents                      <- lapply(parents_children_list,\(pc) pc[["parents"]])
  children                     <- lapply(parents_children_list,\(pc) pc[["children"]])
  block_names                  <- lapply(parents_children_list,\(pc) pc[["names"]])
  block_groups                 <- lapply(parents_children_list,\(pc) pc[["groups"]])


  suppressMessages(simdata_in_list <- lapply(1:length(simdata_list),\(i) {
    simdata_in <- coords_blocking_list[[i]] %>%
      dplyr::select(-na_which) %>%
      dplyr::left_join(simdata_list[[i]])
  }))

  x_list <- lapply(simdata_in_list,\(simdata_in) {
    matrix(0,nrow=nrow(coords),1)
  })


  # remove x_list completely later
  # x_list <- lapply(simdata_in_list,\(simdata_in) {
  #   x <- simdata_in %>%
  #     dplyr::select(dplyr::contains("X_")) %>%
  #     as.matrix()
  #   # colnames(x) <- orig_X_colnames
  #   x[is.na(x)] <- 0 # NAs if added coords due to empty blocks
  #   x
  # })

  # should be common to all images
  indexing_list <- lapply(simdata_in_list,\(simdata_in) {
    blocking <- simdata_in$block %>%
      factor() %>% as.integer()
    indexing <- (1:nrow(simdata_in)-1) %>%
      split(blocking)
  })

  # priors and starting values

  # remove nu dependence later
  start_nu <- 0.5
  matern_fix_twonu <- 1
  # prior and starting values for mcmc

  if(is.null(prior$phi)){
    stop("Need to specify the limits on the Uniform prior for phi via prior$phi.")
  }
  phi_limits <- prior$phi
  if(is.null(starting$phi)){
    start_phi <- mean(phi_limits)
  } else {
    start_phi <- starting$phi
  }

  if(is.null(prior$beta)){
    beta_Vi <- diag(ncol(x_list[[1]])) * 1/100
  } else {
    beta_Vi <- prior$beta
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
    theta_names <- c("phi", "sigmasq")
    npar <- length(theta_names)

    start_theta <- matrix(0, ncol=k, nrow=npar)

    set_unif_bounds <- matrix(0, nrow=npar*k, ncol=2)
    # phi
    set_unif_bounds[seq(1, npar*k, npar),] <- matrix(phi_limits,nrow=1) %x% matrix(1, nrow=k)
    start_theta[1,] <- start_phi

    # sigmasq expansion
    set_unif_bounds[seq(2, npar*k, npar),] <- matrix(c(btmlim, toplim),nrow=1) %x% matrix(1, nrow=k)
    if((q>1)){
      # multivariate without expansion: sigmasq=1 fixed.
      start_theta[2,] <- 1
    } else {
      start_theta[2,] <- btmlim + 1
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
    diag(start_lambda) <- 1
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

  # finally prepare data


  coords <- simdata_in_list[[1]] %>%
    dplyr::select(x,y) %>%
    as.matrix()

  if(verbose > 0){
    cat("Sending to MCMC.\n")
  }

  if(!just_preprocess){
    comp_time <- system.time({
        results <- multi_bipps_mcmc(
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

          # check on removing all of below
          set_unif_bounds,
          beta_Vi,

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

          0, # which_hmc, remove later
          mcmc_adaptive, # adapting

          TRUE, # use_ps, remove later

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

    saved <- listN(y_list, x_list, coords_blocking_list, k,

                   family_id,
      parents, children,
      block_names, block_groups,

      indexing_list,

      set_unif_bounds,
      beta_Vi,


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

  returning <- list(coords = coords,
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


