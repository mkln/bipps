devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

# environment settings
SYSTEM_ENV <- Sys.getenv("SYSTEM_ENV")
if(SYSTEM_ENV == "laptop") {
  sim_idx <- 1
  file_prefix <- "examples/data/"
  n_threads <- 2
  n_samples <- 500
  n_burnin <- 500
  n_thin <- 1
  do_plots <- TRUE
  num_images <- 3
  k <- 2
  q <- 3
} else {
  args <- commandArgs(trailingOnly = TRUE)
  sim_idx <- as.integer(args[1])
  file_prefix <- "/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingpx"
  n_threads <- 4
  n_samples <- 1000
  n_burnin <- 5000
  n_thin <- 5
  do_plots <- FALSE
  num_images <- 40
  k <- 4
  q <- 10
}


# simset1
szs <- c(12,24,48)

# mcmc and model settings
grid <- expand.grid(sz=szs,sim=1:10)
start_idx <- 1
seed_start <- 26

block_size <- 50
starting <- list(phi = 1)
prior <- list(phi = c(0.1,10))
phi_range <- c(1,3)
chains <- 1
sample_theta <- TRUE
mu <- -3

# simulation settings
print(sim_idx)
sim <- grid$sim[sim_idx]
seed <- seed_start + sim
nx <- ny <- grid$sz[sim_idx]
min_nx <- min_ny <- min(grid$sz)
max_nx <- max_ny <- max(grid$sz)

max_n <- max_nx*max_ny
n <- nx*ny
x_max <- 1919
y_max <- 1439
px_x <- x_max / nx
px_y <- y_max / ny
px_area <- px_x * px_y
px_minx <- x_max / min_nx
px_miny <- y_max / min_ny
max_px_area <- px_minx * px_miny
size_factor <- px_area / max_px_area
p <- 1

coords <- expand.grid(x=seq(0,x_max,length.out=max_nx),y=seq(0,y_max,length.out=max_ny))

coords <- coords / max(x_max,y_max)
c_mat <- as.matrix(coords)
d_coords <- as.matrix(dist(c_mat))

save_file <- paste0(file_prefix,"_",sim_idx,".rds")
save_file_ltb <- paste0(file_prefix,"_",sim_idx,"_ltb.rds")

set.seed(seed)
philist <- runif(k,phi_range[1],phi_range[2])

LClist <- 1:k %>% lapply(\(i) t(chol(
  exp(- philist[i] * d_coords)
)))

set.seed(seed)
VV <- lapply(1:num_images,\(j) {
  vlist <- lapply(1:k,\(i) LClist[[i]] %*% rnorm(max_n))

  # factor matrix
  do.call(cbind, vlist)
})

# factor loadings
Lambda <- matrix(0, q, k)
set.seed(seed)
diag(Lambda) <- runif(k, 0.5, 1)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -0.7, 0.7)

Beta <- matrix(rep(mu,p*q),ncol=q)

x_list <- lapply(1:num_images,\(i) {
  matrix(1,nrow = max_n,ncol = p)
})

WW <- lapply(1:num_images,\(i) {
  mat <- VV[[i]] %*% t(Lambda) + x_list[[i]] %*% Beta
  print(range(exp(mat)))
  mat
})

set.seed(seed)
y_agg <- lapply(1:num_images,\(i) {
  yy_list <- lapply(1:k,\(j) {
    mat <- matrix(WW[[i]][,j],nrow = max_ny,ncol = max_nx)
    mat <- matrix(rpois(nrow(mat)*ncol(mat),exp(mat)),nrow=nrow(mat),ncol=ncol(mat))
    r <- raster::raster(mat)
    r_sub <- raster::aggregate(r,fact = max_nx / nx, fun = sum, expand = FALSE)
    vec_sub <- c(raster::as.matrix(r_sub))
  })
  yy <- do.call(cbind,yy_list)
})
agg_coords <- expand.grid(x=seq(0,x_max,length.out=nx),y=seq(0,y_max,length.out=ny))
agg_coords <- agg_coords / max(x_max,y_max)

x_agg <- lapply(1:num_images,\(i) {
  matrix(1,nrow = n,ncol = p)
})

if(do_plots) {
  fields::image.plot(seq(0,x_max,length.out=max_nx),seq(0,y_max,length.out=max_ny),matrix(VV[[1]][,2],nrow = max_ny,ncol = max_nx))
  fields::image.plot(seq(0,x_max,length.out=max_nx),seq(0,y_max,length.out=max_ny),matrix(WW[[1]][,2],nrow = max_ny,ncol = max_nx))
  p1 <- plot_y_list(y_agg,agg_coords)
  print(p1[[1]])
}

# run bipps
out <- lapply(1:chains,\(i) multi_bipps(y_agg,
                                        x_agg,
                                        agg_coords,
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
                                          sample_theta = sample_theta, sample_w = T, sample_lambda = T,
                                          verbose = F, debug = F
                                        ),
                                        just_preprocess = F))
saveRDS(out,save_file)

out_ltb <- lapply(out,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc,beta_mcmc=o$beta_mcmc,waic=o$waic))

saveRDS(out_ltb,save_file_ltb)
