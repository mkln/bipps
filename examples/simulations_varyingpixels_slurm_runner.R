devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
sim_idx <- as.integer(args[1])

# simset1
szs <- c(15,30,45)

# mcmc and model settings
grid <- expand.grid(sz=szs,sim=1:10)
start_idx <- 1
file_prefix <- "/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingpx"
seed_start <- 26

n_samples <- 1000
n_burnin <- 5000
n_thin <- 5
n_threads <- 4
block_size <- 50
starting <- list(phi = 1)
prior <- list(phi = c(0.1,10))
phi_range <- c(1,3)
chains <- 1
do_plots <- FALSE
sample_theta <- TRUE
num_images <- 40
mu <- -0.5
k <- 4
q <- 10

# simulation settings
print(sim_idx)
seed <- seed_start + sim_idx
nx <- ny <- grid$sz[sim_idx]
min_nx <- min_ny <- min(grid$sz)
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

coords <- expand.grid(x=seq(0,x_max,length.out=nx),y=seq(0,y_max,length.out=ny))

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
  wlist <- lapply(1:k,\(i) LClist[[i]] %*% rnorm(n))

  # factor matrix
  do.call(cbind, wlist)
})


# factor loadings
Lambda <- matrix(0, q, k)
set.seed(seed)
diag(Lambda) <- runif(k, 0.5, 1)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -0.7, 0.7)

Beta <- matrix(rep(mu,p*q),ncol=q)

x_list <- lapply(1:num_images,\(i) {
  matrix(1,nrow = n,ncol = p)
})

WW <- lapply(1:num_images,\(i) {
  mat <- VV[[i]] %*% t(Lambda) + x_list[[i]] %*% Beta
  # print(range(exp(mat)))
  mat
})

set.seed(seed)
y_list <- lapply(WW,\(ww) {
  matrix(rpois(nrow(ww)*ncol(ww),exp(ww)*size_factor),nrow=nrow(ww),ncol=ncol(ww))
})

if(do_plots) {
  fields::image.plot(seq(0,x_max,length.out=nx),seq(0,y_max,length.out=ny),matrix(VV[[10]][,2],nrow = ny,ncol = nx))
  fields::image.plot(seq(0,x_max,length.out=nx),seq(0,y_max,length.out=ny),matrix(WW[[10]][,2],nrow = ny,ncol = nx))
  p1 <- plot_y_list(y_list,coords)
  p1[[34]]
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
                                          sample_theta = sample_theta, sample_w = T, sample_lambda = T,
                                          verbose = F, debug = F
                                        ),
                                        just_preprocess = F))
saveRDS(out,save_file)

out_ltb <- lapply(out,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc,beta_mcmc=o$beta_mcmc,waic=o$waic))

saveRDS(out_ltb,save_file_ltb)
