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
  y_list,
  x_list,
  coords,
  k=NULL,
  family = "poisson",
  axis_partition = NULL,
  block_size = 30,
  n_samples = 1000,
  n_burnin = 100,
  n_thin = 1,
  n_threads = 4,
  verbose = 0,
  settings = list(adapting=TRUE,
                    ps=TRUE, saving=TRUE, low_mem=TRUE),
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
){

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

  # data management pt 1
  mcmc_keep <- n_samples
  mcmc_burn <- n_burnin
  mcmc_thin <- n_thin

  which_hmc <- 0

  mcmc_adaptive    <- settings$adapting %>% set_default(TRUE)
  mcmc_verbose     <- debug$verbose %>% set_default(FALSE)
  mcmc_debug       <- debug$debug %>% set_default(FALSE)
  saving <- settings$saving %>% set_default(TRUE)
  low_mem <- settings$low_mem %>% set_default(TRUE)

  debugdag <- debug$dag %>% set_default(1)

  coords %<>% as.matrix()


  # check dimensions of coords and y_list/x_list
  stopifnot("nrow(coords) != nrow(x)"=all(unlist(lapply(x_list,\(x) nrow(x) == nrow(coords)))))
  stopifnot("nrow(coords) != nrow(y)"=all(unlist(lapply(y_list,\(y) nrow(y) == nrow(coords)))))
  stopifnot("ncol(y) not all equal"=all(unlist(lapply(y_list,\(y) ncol(y) == ncol(y_list[[1]])))))


  dd             <- ncol(coords)
  p              <- ncol(x_list[[1]])

  # data management part 0 - reshape/rename
  num_images <- length(y_list)

  y_list <- lapply(y_list,\(y) {
    if(is.null(dim(y))){
      y <- matrix(y, ncol=1)
      orig_y_colnames <<- colnames(y) <- "Y_1"
    } else {
      if(is.null(colnames(y))){
        orig_y_colnames <<- colnames(y) <- paste0('Y_', 1:ncol(y))
      } else {
        orig_y_colnames <<- colnames(y)
        colnames(y)     <- paste0('Y_', 1:ncol(y))
      }
    }
    y
  })


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

  x_list <- lapply(x_list,\(x) {
    if(is.null(colnames(x))){
      orig_X_colnames <<- colnames(x) <- paste0('X_', 1:ncol(x))
    } else {
      orig_X_colnames <<- colnames(x)
      colnames(x)     <- paste0('X_', 1:ncol(x))
    }
    x
  })


  if(is.null(colnames(coords))){
    orig_coords_colnames <- colnames(coords) <- paste0('Var', 1:dd)
  } else {
    orig_coords_colnames <- colnames(coords)
    colnames(coords)     <- paste0('Var', 1:dd)
  }

  q              <- ncol(y_list[[1]])
  k              <- ifelse(is.null(k), q, k)

  # family id
  family <- if(length(family)==1){rep(family, q)} else {family}
  use_ps <- settings$ps %>% set_default(TRUE)
  family_in <- data.frame(family=family)
  available_families <- data.frame(id=c(1,4), family=c("poisson", "negbinomial"))

  suppressMessages(family_id <- family_in %>%
                      left_join(available_families, by=c("family"="family")) %>% pull(id))


  nr             <- nrow(coords)

  ## check this
  if(length(axis_partition) == 1){
    axis_partition <- rep(axis_partition, dd)
  }
  if(is.null(axis_partition)){
    axis_partition <- rep(round((nr/block_size)^(1/dd)), dd)
  }

  # what are we sampling
  sample_w       <- debug$sample_w %>% set_default(TRUE)
  sample_beta    <- debug$sample_beta %>% set_default(TRUE)
  sample_tausq   <- debug$sample_tausq %>% set_default(TRUE)
  sample_theta   <- debug$sample_theta %>% set_default(TRUE)
  sample_lambda  <- debug$sample_lambda %>% set_default(TRUE)


  # data management pt 2

  all_na_which <- lapply(y_list, \(y) apply(y, 1, \(i) ifelse(sum(is.na(i))==q,NA,1)))

  fixed_thresholds <- 1:dd %>% lapply(function(i) kthresholdscp(coords[,i], axis_partition[i]))

  # Domain partitioning and gibbs groups

  system.time(coords_blocking_list <- lapply(all_na_which,\(na_which_i) {
    coords_blocking <- coords %>%
                  as.matrix() %>%
                  tessellation_axis_parallel_fix(fixed_thresholds, 1) %>%
                  dplyr::mutate(na_which = na_which_i)
    coords_blocking
  }))

  # this depends on there being no all NA blocks, otherwise the coords_blocking_list
  # will not be the same across images. Change this?
  # DAG
  system.time(
    if(dd < 4){
      pc <- mesh_graph_build(coords_blocking_list[[1]], axis_partition, FALSE, n_threads, debugdag)
    } else {
      stop("Input dimension is too high?!")
    }
  )



  parents                      <- pc[["parents"]]
  children                     <- pc[["children"]]
  block_names                  <- pc[["names"]]
  block_groups                 <- pc[["groups"]]

  if(indpart){
    parents %<>% lapply(function(x) lapply(x,\(x_i) numeric(0)))
    children %<>% lapply(function(x) numeric(0))
    block_groups %<>% rep(0, length(block_groups))
  }

  blocking <- coords_blocking_list[[1]]$block %>%
      factor() %>% as.integer()
  indexing <- (1:nrow(coords)-1) %>%
      split(blocking)

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
      start_v <- lapply(1:num_images,\(i) matrix(0, nrow = nrow(coords), ncol = k))
    } else {
      # this is used to restart MCMC
      # assumes the ordering and the sizing is correct,
      # so no change is necessary and will be input directly to mcmc
      start_v <- starting$v
    }
  }



  coords_renamer <- colnames(coords)
  names(coords_renamer) <- orig_coords_colnames


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
          indexing,

          set_unif_bounds,
          beta_Vi,


          # sigmasq_ab,
          # tausq_ab,

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

    saved <- listN(y_list, x_list, coords_blocking_list, k,

                   family_id,
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

