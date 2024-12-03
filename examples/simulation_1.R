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

nx <- 1000
ny <- 1000
num_images <- 3
# theta
theta <- 0.7 # scaled version
inv_theta <- 1 / theta * max(x_max,y_max)
sigmasq <- 1
x_max <- 1919
y_max <- 1439
scaling <- 20
mu <- -9
k <- 3
q <- 10
p <- 1

# generate point process
# pp <- rLGCP(model="exponential",param=list(var=sigmasq,scale=theta),win=owin(c(0,1919),c(0,1439)))
# pp
# plot(pp)

grid<- list( x= seq( 0,x_max*scaling,,nx), y= seq(0,y_max*scaling,,ny))
obj<- circulantEmbeddingSetup( grid, Covariance="Exponential", theta=inv_theta)

# subset
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

    # dimx <- gridx[2] - gridx[1]
    # dimy <- gridy[2] - gridy[1]
    #
    #
    # # px_ar <- dimx * dimy
    #
    # exp_sim <- exp(subW)
  })
})

out1 <- readRDS("out1_small_sample.rds")[[1]]

thidx <- seq(1,800,by=40)
beta0 <- out1$plus_icept_mcmc[,thidx]
rowMeans(beta0)

v_mcmc <- out1$v_mcmc
v_mcmc <- simplify2array(v_mcmc)
range(apply(v_mcmc,1:2,mean))

lambda <- get_rvars(list(out1),"lambda",thin=40)
theta <- get_rvars(list(out1),"theta",thin=40)


image.plot( gridx,gridy, exp(V[[1]][[1]]))


VV <- lapply(1:num_images,\(i) {
  lapply(V[[i]],\(w) c(w)) %>%
    do.call(cbind,.)
})


# factor loadings
Lambda <- matrix(0, q, k)
diag(Lambda) <- runif(k, 0.5, 2)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -1, 1)

# Lambda <- diag(q)
# Lambda[2,2] <- .5

# nuggets
# tau.sq <- rep(.01, q)
# TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
# # measurement errors
# # double check these
# EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

XX <- lapply(1:num_images,\(j) {
  1:p %>% lapply(function(i) rnorm(nrow(VV[[i]]), 0, .1)) %>% do.call(cbind, .)
})
Beta <- matrix(rnorm(p*q), ncol=q) * 0

WW <- lapply(1:num_images,\(i) VV[[i]] %*% t(Lambda))

sz_x <- length(gridx)
sz_y <- length(gridy)
W <- lapply(1:num_images,\(i) {
  lapply(1:q,\(j) {
    mat <- matrix(WW[[i]][,j],nrow=sz_y,ncol=sz_x)
    print(range(mat))
    mat
  })
})

lin_pred <- lapply(1:num_images,\(i) {
  lapply(1:q,\(j) {
    mat <- W[[i]][[j]] + mu #+ matrix(XX[[i]][,j] %*% Beta,nrow=sz_y,ncol=sz_x)
    print(range(mat))
    mat
  })
})

image.plot( gridx,gridy, exp(V[[2]][[1]]))
image.plot( gridx,gridy, exp(W[[2]][[1]]))
image.plot( gridx,gridy, exp(lin_pred[[2]][[1]]))


window <- owin(c(0,gridx[length(gridx)]),c(0,gridy[length(gridy)]))

imgs <- lapply(1:num_images,\(i) {
  lapply(1:q,\(j) {
    img <- as.im(t(exp(lin_pred[[i]][[j]])),W=window)
    print(mean(img)*x_max*y_max)
    img
  })
})

img <- imgs[[1]][[1]]
mean(img)*x_max*y_max
#
pats <- lapply(1:num_images,\(i) {
  pp <- rmpoispp(imgs[[i]])
})

pp_df <- lapply(1:num_images,\(i) {
  pats[[i]] %>%
    as.data.frame() %>%
    rename(type = marks) %>%
    mutate(spot = paste0("spot_",i))
}) %>%
  bind_rows()

pp_df %>%
  ggplot(aes(x,y,color=type)) +
  geom_point() +
  theme_classic() +
  facet_wrap(~spot,ncol=1)

out <- pixellate_grid(pp_df$x,pp_df$y,pp_df$type,pp_df$spot,28,21)

y_list <- out$y_list
coords <- out$coords
x_list <- lapply(y_list,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

p1 <- plot_y_list(y_list,coords)
p1

n_samples <- 20
n_burnin <- 1000
n_thin <- 1
n_threads <- 4
block_size <- 50
starting <- list(phi = 5)
# D <- sqrt(max_dim[1]^2 + max_dim[2]^2)
# a <- -log(0.05) * 3 / D
prior <- list(phi = c(0.1,10))

chains <- 1

out1 <- lapply(1:chains,\(i) multi_bipps(y_list,
                                         x_list,
                                         coords,
                                         k = k,
                                         family = "poisson",
                                         block_size = block_size,
                                         n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                                         n_threads = n_threads,
                                         starting = starting,
                                         prior = prior,
                                         settings = list(adapting = T, saving = T, ps = T,low_mem=F),
                                         verbose = 10,
                                         debug = list(
                                           sample_beta = T, sample_tausq = F,
                                           sample_theta = T, sample_w = T, sample_lambda = T,
                                           verbose = F, debug = F
                                         ),
                                         just_preprocess = F))

out <- readRDS("out_sim1.rds")

lambda <- get_rvars(out,"lambda",thin=n_thin)
lambda
theta <- get_rvars(out,"theta",thin=n_thin)
theta
mcmc_trace(as_draws_df(theta[1,]))

