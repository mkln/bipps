## SIM STUDY 1
# 1. 2 outcomes, 2 latent factors
# 2. varying number of images and grid cells
# 3. fix MCMC samples at 5000 burn in, 1000 samples
# 4. sample everything except tau^2
# 5.

devtools::load_all()
library(tidyverse)

grid <- expand.grid(
  num_images=c(2,5,10,15,20,25,30,35,40,45,50),
  num_grids=c(20,30,40),
  rep=seq(1,5)
)

fit_multi_bipps <- function(num_images,num_grids,rep) {
  set.seed(2020+rep)
  SS <- num_grids
  dd <- 2 # spatial dimension
  n <- SS^2 # number of locations
  q <- 2 # number of outcomes
  k <- 2 # number of spatial factors used to make the outcomes
  p <- 1 # number of covariates

  xlocs <- seq(0, 1, length.out=SS)
  coords <- expand.grid(list(xlocs, xlocs)) %>%
    as.data.frame()

  clist <- 1:q %>% lapply(function(i) coords %>%
                            mutate(mv_id=i) %>%
                            as.matrix())

  philist <- c(1, 1) # spatial decay for each factor

  # cholesky decomp of covariance matrix
  LClist <- 1:k %>% lapply(function(i) t(chol(
    #exp(- philist[i] * as.matrix(dist(clist[[i]])))))) #^2 + diag(nrow(clist[[i]]))*1e-5))))
    bipps:::Cov_matern(clist[[i]], clist[[i]], 1, philist[i], 0.5, 0, T, 10))))

  # generating the factors
  WW <- lapply(1:num_images,\(j) {
    wlist <- 1:k %>% lapply(function(i) LClist[[i]] %*% rnorm(n))

    # factor matrix
    do.call(cbind, wlist)
  })


  # factor loadings
  Lambda <- matrix(0, q, ncol(WW[[1]]))
  diag(Lambda) <- runif(k, 1, 2)
  Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -1, 1)
  Lambda[2,2] <- .5

  # nuggets
  # tau.sq <- rep(.01, q)
  # TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
  # # measurement errors
  # # double check these
  # EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

  XX <- lapply(1:num_images,\(j) {
    1:p %>% lapply(function(i) rnorm(n, 0, .1)) %>% do.call(cbind, .)
  })
  Beta <- matrix(rnorm(p*q), ncol=q) * 0

  # outcome matrix, fully observed
  linear_predictor <- lapply(1:num_images,\(i) XX[[i]] %*% Beta + WW[[i]] %*% t(Lambda))
  YY <- lapply(1:num_images,\(i) {
    y <- matrix(0, ncol=q, nrow=nrow(linear_predictor[[i]]))
    y <- lapply(1:q,\(j) {
      rpois(n, exp(linear_predictor[[i]][,j]))
    })
    y <- do.call(cbind,y)
  })


  mcmc_keep <- 1000
  mcmc_burn <- 5000
  mcmc_thin <- 1

  # Lambda_start <- matrix(0, q, ncol(WW[[1]]))\
  Lambda_start <- Lambda
  Beta_start <- Beta

  mesh_total_time <- system.time({
    meshout <- multi_bipps(YY, family="poisson", XX, coords, k = 2,
                           block_size=25,
                           n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin,
                           n_threads = 1,
                           starting=list(lambda = Lambda_start, beta=Beta_start, phi=1),
                           prior = list(btmlim= .01, toplim=1e3, phi=c(.1, 20), nu=c(.5, .5)),
                           settings = list(adapting=T, saving=F, ps=T, hmc=0),
                           verbose=10,
                           debug=list(sample_beta=T, sample_tausq=F,
                                      sample_theta=T, sample_w=T, sample_lambda=T,
                                      verbose=F, debug=F)
    )})

  list(result=meshout,WW=WW,YY=YY,Lambda=Lambda,num_images=num_images,num_grids=num_grids,rep=rep)
}

out <- fit_multi_bipps(4,20,3)
out$result$lambda_mcmc %>% apply(1:2, mean)
out$result$theta_mcmc %>% apply(1:2, mean)
out$Lambda
(lower <- out$result$theta_mcmc %>%
    apply(., 1:2, \(x) quantile(x, probs = c(0.025))))

(upper <- out$result$theta_mcmc %>%
    apply(., 1:2, \(x) quantile(x, probs = c(0.975))))

# need a decent estimate of lambda, I think, in order to correctly estimate theta and lambda.
# How can we get a good estimate of lambda to start out with?
# Probably the easiest is just to fit the model at low-res quickly and then start with that estimate.

rslurm::slurm_apply(fit_multi_bipps,grid,jobname="bipps_sim_study_1",nodes=dim(grid)[1],cpus_per_node=1,submit=FALSE,upload="joelne",
                    slurm_options = list(mem="8G",time="20:00"))
