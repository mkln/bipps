rm(list=ls())
devtools::load_all()
library(magrittr)
library(dplyr)
library(ggplot2)

set.seed(2020)

SS <- 20 # coord values for jth dimension
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

philist <- c(15, 15) # spatial decay for each factor

# cholesky decomp of covariance matrix
LClist <- 1:k %>% lapply(function(i) t(chol(
  #exp(- philist[i] * as.matrix(dist(clist[[i]])))))) #^2 + diag(nrow(clist[[i]]))*1e-5))))
  bipps:::Cov_matern(clist[[i]], clist[[i]], 1, philist[i], 0.5, 0, T, 10))))

# generating the factors
wlist <- 1:k %>% lapply(function(i) LClist[[i]] %*% rnorm(n))

# factor matrix
WW <- do.call(cbind, wlist) * 1e-3

# factor loadings
Lambda <- matrix(0, q,ncol(WW))
diag(Lambda) <- runif(k, 1, 2)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -1, 1)
Lambda[2,2] <- .5

# nuggets
tau.sq <- rep(.01, q)
TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
# measurement errors
EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

XX <- 1:p %>% lapply(function(i) rnorm(n, 0, .1)) %>% do.call(cbind, .)
Beta <- matrix(rnorm(p*q), ncol=q) * 0

# outcome matrix, fully observed
linear_predictor <- XX %*% Beta + WW %*% t(Lambda)
YY_full <- matrix(0, ncol=q, nrow=nrow(linear_predictor))

YY_full[,1] <- rpois(n, exp(linear_predictor[,1]))
YY_full[,2] <- rpois(n, exp(linear_predictor[,2]))

# .. introduce some NA values in the outcomes
YY <- YY_full

#YY[sample(1:n, n/5, replace=FALSE), 1] <- NA
#YY[sample(1:n, n/5, replace=FALSE), 2] <- NA


simdata <- coords %>%
  cbind(data.frame(Outcome_full=YY_full,
                   # Outcome_obs = YY,
                   linear_pred = linear_predictor))

simdata %>%
  tidyr::gather(Outcome, Value, -all_of(colnames(coords))) %>%
  ggplot(aes(Var1, Var2, fill=Value)) +
  geom_raster() +
  facet_wrap(Outcome ~., ncol=2, scales="free") +
  scale_fill_viridis_c()

mcmc_keep <- 1000
mcmc_burn <- 5000
mcmc_thin <- 1


set.seed(1)
mesh_total_time <- system.time({
  meshout <- bipps(YY, family="poisson", XX, coords, k = 2,
                      block_size=25,
                      n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin,
                      n_threads = 16,
                      starting=list(lambda = Lambda, beta=Beta, phi=1),
                      prior = list(btmlim= .01, toplim=1e3, phi=c(.1, 20), nu=c(.5, .5)),
                      settings = list(adapting=T, saving=F, ps=T, hmc=0),
                      verbose=10,
                      debug=list(sample_beta=T, sample_tausq=F,
                                 sample_theta=T, sample_w=T, sample_lambda=T,
                                 verbose=F, debug=F)
  )})


plot_cube <- function(cube_mcmc, q, k, name="Parameter"){
  par(mar=c(2.5,2,1,1), mfrow=c(q,k))
  for(i in 1:q){
    for(j in 1:k){
      cube_mcmc[i, j,] %>% plot(type='l', main="{name} {i}, {j}" %>% glue::glue())
    }
  }
}

# chain plots
plot_cube(meshout$theta_mcmc, 1, k, "theta")
plot_cube(meshout$lambda_mcmc, q, k, "Lambda")
plot_cube(meshout$beta_mcmc, p, q, "Beta")

# posterior means
meshout$tausq_mcmc %>% apply(1, mean)
meshout$lambda_mcmc %>% apply(1:2, mean)
meshout$beta_mcmc %>% apply(1:2, mean)
meshout$theta_mcmc %>% apply(1:2, mean)

# process means
wmesh <- data.frame(meshout$w_mcmc %>% summary_list_mean())
colnames(wmesh) <- paste0("wmesh_", 1:k)
# predictions
ymesh <- data.frame(meshout$yhat_mcmc %>% summary_list_mean())
colnames(ymesh) <- paste0("ymesh_", 1:q)

mesh_df <-
  meshout$coordsdata %>%
  cbind(ymesh,wmesh)
results <- simdata %>% left_join(mesh_df)

# prediction rmse, out of sample
results %>% filter(!complete.cases(Outcome_obs.1)) %>%
  with((Outcome_full.1 - ymesh_1)^2) %>% mean() %>% sqrt()
results %>% filter(!complete.cases(Outcome_obs.2)) %>%
  with((Outcome_full.2 - ymesh_2)^2) %>% mean() %>% sqrt()

(postmeans <- results %>% dplyr::select(Var1, Var2, ymesh_1) %>%
    tidyr::gather(Variable, Value, -Var1, -Var2) %>%
    ggplot(aes(Var1, Var2, fill=Value)) +
    geom_raster() +
    facet_wrap(Variable ~ ., ncol= 2) +
    scale_fill_viridis_c())


