# spatially stationary W
# (also try with nonstationary)
# vary -> number of outcomes, spatial factors and images. lambda and theta.

devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)
library(fields)

set.seed(2020)

# mcmc settings
n_samples <- 500
n_burnin <- 5000
n_thin <- 10
n_threads <- 16
block_size <- 50
starting <- list(phi = 5)
prior <- list(phi = c(0.1,10))
chains <- 2
do_plots <- FALSE
save_file <- "out_sim2.rds"
save_file_lt <- "out_sim2_lt.rds"

# simulation settings
nx <- 30
ny <- 30
num_images <- 30
# theta
x_max <- 1919
y_max <- 1439
Theta <- 0.7 # scaled version
inv_theta <- 1 / Theta * max(x_max,y_max)
sigmasq <- 1

scaling <- 20
mu <- -1
k <- 3
q <- 5
p <- 1

grid<- list( x= seq( 0,x_max*scaling,,nx*scaling), y= seq(0,y_max*scaling,,ny*scaling))
obj<- circulantEmbeddingSetup( grid, Covariance="Exponential", theta=inv_theta)

# subset V
idx_x <- which(grid[[1]] < x_max)
idx_y <- which(grid[[2]] < y_max)

gridx <- grid[[1]][idx_x]
gridy <- grid[[2]][idx_y]

set.seed(2025)
V <- lapply(1:num_images,\(i) {
  lapply(1:k,\(i) {
    sim<-  sqrt(sigmasq)*circulantEmbedding( obj)
    subW <- sim[idx_y,idx_x]

    print(range(subW))
    subW
  })
})


VV <- lapply(1:num_images,\(i) {
  lapply(V[[i]],\(w) c(w)) %>%
    do.call(cbind,.)
})


# factor loadings
set.seed(2025)
Lambda <- matrix(0, q, k)
diag(Lambda) <- runif(k, 0.5, 1)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -0.7, 0.7)

# nuggets
# tau.sq <- rep(.01, q)
# TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
# # measurement errors
# # double check these
# EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

# XX <- lapply(1:num_images,\(j) {
#   1:p %>% lapply(function(i) rnorm(nrow(VV[[i]]), 0, .1)) %>% do.call(cbind, .)
# })
# Beta <- matrix(rnorm(p*q), ncol=q) * 0

WW <- lapply(1:num_images,\(i) {
  mat <- VV[[i]] %*% t(Lambda) + mu
  print(range(exp(mat)))
  mat
})
# sz_x <- length(gridx)
# sz_y <- length(gridy)
# W <- lapply(1:num_images,\(i) {
#   lapply(1:q,\(j) {
#     mat <- matrix(WW[[i]][,j],nrow=sz_y,ncol=sz_x)
#     print(range(mat))
#     mat
#   })
# })
#
# # make linear predictor
# lin_pred <- lapply(1:num_images,\(i) {
#   lapply(1:q,\(j) {
#     mat <- W[[i]][[j]] + mu #+ matrix(XX[[i]][,j] %*% Beta,nrow=sz_y,ncol=sz_x)
#     print(range(mat))
#     mat
#   })
# })

y_list <- lapply(WW,\(ww) {
  matrix(rpois(nrow(ww)*ncol(ww),exp(ww)),nrow=nrow(ww),ncol=ncol(ww))
})

coords <- expand.grid(x=gridx,y=gridy)

coords_scaled <- coords / max(x_max,y_max)


# make point patterns
# window <- owin(c(0,gridx[length(gridx)]),c(0,gridy[length(gridy)]))
#
# imgs <- lapply(1:num_images,\(i) {
#   lapply(1:q,\(j) {
#     img <- as.im(t(exp(lin_pred[[i]][[j]])),W=window)
#     print(mean(img)*x_max*y_max)
#     img
#   })
# })

# img <- imgs[[1]][[1]]
# mean(img)*x_max*y_max
#
# pats <- lapply(1:num_images,\(i) {
#   pp <- rmpoispp(imgs[[i]])
# })
#
# pp_df <- lapply(1:num_images,\(i) {
#   pats[[i]] %>%
#     as.data.frame() %>%
#     rename(type = marks) %>%
#     mutate(spot = paste0("spot_",i))
# }) %>%
#   bind_rows()

# pp_df %>%
#   ggplot(aes(x,y,color=type)) +
#   geom_point() +
#   theme_classic() +
#   facet_wrap(~spot,ncol=1)

# pixellate the point patterns
# out <- pixellate_grid(pp_df$x,pp_df$y,pp_df$type,pp_df$spot,28,21)
#
# y_list <- out$y_list
# coords <- out$coords
x_list <- lapply(y_list,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

if(do_plots) {
  p1 <- plot_y_list(y_list,coords_scaled)
  p1
}

# run bipps
out <- lapply(1:chains,\(i) multi_bipps(y_list,
                                         x_list,
                                         coords,
                                         k = k,
                                         family = "poisson",
                                         block_size = block_size,
                                         n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                                         n_threads = n_threads,
                                         starting = starting,
                                         prior = prior,
                                         settings = list(adapting = T, saving = T, ps = T,low_mem=T),
                                         verbose = 10,
                                         debug = list(
                                           sample_beta = T, sample_tausq = F,
                                           sample_theta = T, sample_w = T, sample_lambda = T,
                                           verbose = F, debug = F
                                         ),
                                         just_preprocess = F))
saveRDS(out,save_file)

out_lt <- lapply(out,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc))

saveRDS(out_lt,save_file_lt)
