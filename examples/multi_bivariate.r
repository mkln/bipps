#renv::activate(".")
rm(list=ls())
devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(2020)

num_images <- 5

SS <- 30 # coord values for jth dimension
dd <- 2 # spatial dimension
n <- SS^2 # number of locations
q <- 4 # number of outcomes
k <- 2 # number of spatial factors used to make the outcomes
p <- 1 # number of covariates
theta <- 0.7

xlocs <- seq(0, 1, length.out=SS)
coords <- expand.grid(list(xlocs, xlocs)) %>%
  as.data.frame()

clist <- 1:q %>% lapply(function(i) coords %>%
                          mutate(mv_id=i) %>%
                          as.matrix())

philist <- rep(theta,k) # spatial decay for each factor

# cholesky decomp of covariance matrix
LClist <- 1:k %>% lapply(function(i) t(chol(
  exp(- philist[i] * as.matrix(dist(clist[[i]])))
  ))) #^2 + diag(nrow(clist[[i]]))*1e-5))))
  # bipps:::Cov_matern(clist[[i]], clist[[i]], 1, philist[i], 1.5, 0, T, 10))))

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
# Lambda[2,2] <- .5

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
YY_full <- lapply(1:num_images,\(i) {
  y <- matrix(0, ncol=q, nrow=nrow(linear_predictor[[i]]))
  y <- lapply(1:q,\(j) {
    rpois(n, exp(linear_predictor[[i]][,j]))
  })
  y <- do.call(cbind,y)
})

YY <- YY_full


lapply(1:length(YY),\(i) {
  y <- YY[[i]]
  as_tibble(y) %>%
    mutate(image = i) %>%
    bind_cols(coords)
}) %>%
  bind_rows() %>%
  rename(X = Var1, Y = Var2) %>%
  mutate( image = factor(image)) -> simdata

plotting_data <- TRUE
if(plotting_data){
  simdata %>%
    pivot_longer(c(V1,V2,V3,V4)) %>%
    filter(image %in% c(1:10)) %>%
    mutate(name=factor(name)) %>%
    ggplot(aes(X,Y,fill=log(value+1))) +
    geom_raster() +
    facet_grid(image~name)
}

mcmc_keep <- 1000
mcmc_burn <- 500

mcmc_thin <- 1

na_sz <- 20
x_na <- sample(1:SS^dd,na_sz,replace = TRUE)
y_na <- sample(1:q,na_sz,replace = TRUE)
YY <- lapply(YY,\(yy) {
  yy[cbind(x_na,y_na)] <- NA
  yy
})


#devtools::load_all()
set.seed(1)
mesh_total_time <- system.time({
  meshout <- multi_bipps(YY, family="poisson", XX, coords, k = k,
                      block_size=25,
                      n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin,
                      n_threads = 2,
                      starting=list(lambda = Lambda, beta=Beta, phi=15),
                      prior = list(btmlim= .01, toplim=1e3, phi=c(1, 40), nu=c(1.5, 1.5)),
                      settings = list(adapting=T, saving=F, ps=F, hmc=0),
                      verbose=10,
                      debug=list(sample_beta=T, sample_tausq=F,
                                 sample_theta=T, sample_w=T, sample_lambda=T,
                                 verbose=F, debug=F)
  )})


# plot_cube <- function(cube_mcmc, q, k, name="Parameter"){
#   par(mar=c(2.5,2,1,1), mfrow=c(q,k))
#   for(i in 1:q){
#     for(j in 1:k){
#       cube_mcmc[i, j,] %>% plot(type='l', main="{name} {i}, {j}" %>% glue::glue())
#     }
#   }
# }
#
# # chain plots
# plot_cube(meshout$theta_mcmc, 1, k, "theta")
# plot_cube(meshout$lambda_mcmc, q, k, "Lambda")
# plot_cube(meshout$beta_mcmc, p, q, "Beta")
#
# # posterior means
# meshout$tausq_mcmc %>% apply(1, mean)
# meshout$lambda_mcmc %>% apply(1:2, mean)
# meshout$beta_mcmc %>% apply(1:2, mean)
# meshout$theta_mcmc %>% apply(1:2, mean)
#
# # process means
# wmesh <- data.frame(meshout$w_mcmc %>% summary_list_mean())
# colnames(wmesh) <- paste0("wmesh_", 1:k)
# # predictions
# ymesh <- data.frame(meshout$yhat_mcmc %>% summary_list_mean())
# colnames(ymesh) <- paste0("ymesh_", 1:q)
#
# mesh_df <-
#   meshout$coordsdata %>%
#   cbind(ymesh) %>%
#   # cbind(wmesh) %>%
#   cbind(YY)
#
# results <- simdata %>% left_join(mesh_df)
#
# # prediction rmse, out of sample
# results %>% filter(!complete.cases(Outcome_obs.1)) %>%
#   with((Outcome_full.1 - ymesh_1)^2) %>% mean() %>% sqrt()
# results %>% filter(!complete.cases(Outcome_obs.2)) %>%
#   with((Outcome_full.2 - ymesh_2)^2) %>% mean() %>% sqrt()
#
# (postmeans <- results %>% dplyr::select(Var1, Var2, ymesh_1) %>%
#     tidyr::gather(Variable, Value, -Var1, -Var2) %>%
#     ggplot(aes(Var1, Var2, fill=Value)) +
#     geom_raster() +
#     facet_wrap(Variable ~ ., ncol= 2) +
#     scale_fill_viridis_c())
#
# res <- mesh_df %>%
#   pivot_longer(-c(x,y),names_pattern = "(pred|obs)_(.*)$",names_to=c(".value","type"))
#
#
# res %>%
#   pivot_longer(c(pred,obs)) %>%
#   ggplot(aes(x,y, fill=value)) +
#   geom_tile() +
#   facet_wrap(name ~ type, ncol= q) +
#   scico::scale_fill_scico(palette = "batlow")
#
# res %>%
#   group_by(type) %>%
#   summarise(rmse=sqrt(mean((pred-obs)^2)))
#
# (lower <- meshout$lambda_mcmc %>%
#     apply(.,1:2,\(x) quantile(x,probs=c(0.025))))
#
# (upper <- meshout$lambda_mcmc %>%
#     apply(.,1:2,\(x) quantile(x,probs=c(0.975))))
