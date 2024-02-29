rm(list=ls())

library(magrittr)
library(dplyr)
library(ggplot2)

set.seed(2020)

q <- 5 # number of outcomes
k <- 3 # number of spatial factors used to make the outcomes
p <- 1 # number of covariates

# coordinate system is common for all subjects
# as well as covariance hyperparameters for all spatial processes

SS <- 15 # coord values for jth dimension 
dd <- 2 # spatial dimension
n <- SS^2 # number of locations
xlocs <- seq(0, 1, length.out=SS)
coords <- expand.grid(list(xlocs, xlocs)) %>% as.matrix()

philist <- rep(1, k) # spatial decay for each factor

# cholesky decomp of covariance matrix
LClist <- 1:k %>% lapply(function(i) t(chol(
  meshed:::Cov_matern(coords, coords, 1, philist[i], 0.5, 0, T, 10))))


# factor loadings, common for all subjects
Lambda <- matrix(0, q, k)
diag(Lambda) <- runif(k, 1, 2)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -1, 1)
Lambda[2,2] <- .5

# covariate effects for each outcome 
# !! check that Y values are reasonable because exp(XB) may be very large
Beta <- matrix(rnorm(p*q), ncol=q) 

gen_subj_data <- function(subj, coords, Beta, Lambda, LClist){
  
  q <- ncol(Beta)
  n <- nrow(coords)
  p <- nrow(Beta)
  k <- ncol(Lambda)
  
  # generating the factors
  wlist <- 1:k %>% lapply(function(i) LClist[[i]] %*% rnorm(n))
  
  # factor matrix
  WW <- do.call(cbind, wlist)
  
  XX <- 1:p %>% lapply(function(i) rnorm(n, 0, .1)) %>% do.call(cbind, .)
  
  
  # outcome matrix, fully observed
  LambdaW <- tcrossprod(WW, Lambda)
  linear_predictor <- XX %*% Beta + LambdaW
  YY_full <- matrix(0, ncol=q, nrow=nrow(linear_predictor))
  
  for(j in 1:q){
    YY_full[,j] <- rpois(n, exp(linear_predictor[,j]))
  }
  
  return(cbind(data.frame(coords), data.frame(subject=subj,
              X = XX,
              W = WW,
              LambdaW = LambdaW,
              Y = YY_full)))
}


nsub <- 10
subject_data <- list()
for(i in 1:nsub){
  subject_data[[i]] <- gen_subj_data(i, coords, Beta, Lambda, LClist)
}

df_all_subj <- bind_rows(subject_data)

# for each subject, plot latent intensity for the first outcome (i.e., cell type)
df_all_subj %>%
  ggplot(aes(Var1, Var2, fill=LambdaW.1)) + 
  geom_raster() + 
  facet_wrap(subject ~., ncol=2, scales="free") +
  scale_fill_viridis_c()
