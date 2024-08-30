get_test_data_multi <- function(n_jobs=1,num_images=3) {

  SS <- 16 # coord values for jth dimension
  dd <- 2 # spatial dimension
  n <- SS^2 # number of locations
  q <- 4 # number of outcomes
  k <- 2 # number of spatial factors used to make the outcomes
  p <- 1 # number of covariates

  xlocs <- seq(0, 1, length.out=SS)
  coords <- expand.grid(list(xlocs, xlocs)) %>%
    as.data.frame()

  clist <- 1:q %>% lapply(function(i) coords %>%
                            mutate(mv_id=i) %>%
                            as.matrix())

  philist <- c(15, 15) # spatial decay for each factor

  # cholesky decomp of covariance matrix
  LClist <- 1:k %>% lapply(function(i) t(chol(
    #exp(- philist[i] * as.matrix(dist(clist[[i]])))))) #^2 + diag(nrow(clist[[i]]))*1e-5))))
    bipps:::Cov_matern(clist[[i]], clist[[i]], 1, philist[i], 1.5, 0, T, 10))))

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

  mcmc_keep <- 50
  mcmc_burn <- 10

  mcmc_thin <- 1

  meshout <- multi_bipps(YY, family="poisson", XX, coords, k = 2,
                         block_size=25,
                         n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin,
                         n_threads = n_jobs,
                         starting=list(lambda = Lambda, beta=Beta, phi=15),
                         prior = list(btmlim= .01, toplim=1e3, phi=c(1, 40), nu=c(1.5, 1.5)),
                         settings = list(adapting=T, saving=F, ps=F, hmc=0),
                         verbose=0,
                         debug=list(sample_beta=T, sample_tausq=F,
                                    sample_theta=T, sample_w=T, sample_lambda=T,
                                    verbose=F, debug=F)
  )

  idx <- seq(mcmc_burn,mcmc_keep,5)

  test_data <- list(beta_test=meshout$beta_mcmc[,,idx],
                    lambda_test=meshout$lambda_mcmc[,,idx],
                    theta_test=meshout$theta_mcmc[,,idx],
                    w_test=meshout$w_mcmc[idx])
}
