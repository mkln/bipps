devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(bayesplot)
library(ggdist)
library(patchwork)
library(latex2exp)

theme_set(theme_bw(base_size=11, base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

fsave <- \(fname) {
  ggsave(fname,dpi=300, height=5, width=8, units="in")
}
figures_folder <- "examples/data/figures/varyingpixels/"

# simset1
szs <- c(12,24,48)

# mcmc and model settings
grid <- expand.grid(sz=szs,sim=1:10)
start_idx <- 1
file_prefix <- "examples/data/out_simset_varyingpx"
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
mu <- -3
k <- 4
q <- 10
x_max <- 1919
y_max <- 1439
p <- 1


unique_combinations_with_self <- function(data) {
  # Generate all pairs including self-pairings
  all_pairs <- expand.grid(data, data)

  # Sort the pairs and remove duplicates
  all_pairs <- t(apply(all_pairs, 1, sort))
  unique_pairs <- unique(all_pairs)

  return(unique_pairs)
}

# uncomment this to run the first time, save, and then comment out again
# lapply(start_idx:nrow(grid),\(sim_idx) {
#   print(sim_idx)
#   sim <- grid$sim[sim_idx]
#   seed <- seed_start + sim
#   nx <- ny <- grid$sz[sim_idx]
#   min_nx <- min_ny <- min(grid$sz)
#   max_nx <- max_ny <- max(grid$sz)
#
#   max_n <- max_nx*max_ny
#   n <- nx*ny
#   x_max <- 1919
#   y_max <- 1439
#   px_x <- x_max / nx
#   px_y <- y_max / ny
#   px_area <- px_x * px_y
#   px_minx <- x_max / min_nx
#   px_miny <- y_max / min_ny
#   max_px_area <- px_minx * px_miny
#   size_factor <- px_area / max_px_area
#   p <- 1
#
#   coords <- expand.grid(x=seq(0,x_max,length.out=max_nx),y=seq(0,y_max,length.out=max_ny))
#
#   coords <- coords / max(x_max,y_max)
#   c_mat <- as.matrix(coords)
#   d_coords <- as.matrix(dist(c_mat))
#
#   save_file <- paste0(file_prefix,"_",sim_idx,".rds")
#   save_file_ltb <- paste0(file_prefix,"_",sim_idx,"_ltb.rds")
#
#   set.seed(seed)
#   philist <- runif(k,phi_range[1],phi_range[2])
#
#   LClist <- 1:k %>% lapply(\(i) t(chol(
#     exp(- philist[i] * d_coords)
#   )))
#
#   set.seed(seed)
#   VV <- lapply(1:num_images,\(j) {
#     vlist <- lapply(1:k,\(i) LClist[[i]] %*% rnorm(max_n))
#
#     # factor matrix
#     do.call(cbind, vlist)
#   })
#
#   # factor loadings
#   Lambda <- matrix(0, q, k)
#   set.seed(seed)
#   diag(Lambda) <- runif(k, 0.5, 1)
#   Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -0.7, 0.7)
#
#   Beta <- matrix(rep(mu,p*q),ncol=q)
#
#   x_list <- lapply(1:num_images,\(i) {
#     matrix(1,nrow = max_n,ncol = p)
#   })
#
#   WW <- lapply(1:num_images,\(i) {
#     mat <- VV[[i]] %*% t(Lambda) + x_list[[i]] %*% Beta
#     print(range(exp(mat)))
#     mat
#   })
#
#   set.seed(seed)
#   y_list <- lapply(1:num_images,\(i) {
#     yy_list <- lapply(1:q,\(j) {
#       mat <- matrix(WW[[i]][,j],nrow = max_ny,ncol = max_nx)
#       mat <- matrix(rpois(nrow(mat)*ncol(mat),exp(mat)),nrow=nrow(mat),ncol=ncol(mat))
#       r <- raster::raster(mat)
#       r_sub <- raster::aggregate(r,fact = max_nx / nx, fun = sum, expand = FALSE)
#       vec_sub <- c(raster::as.matrix(r_sub))
#     })
#     yy <- do.call(cbind,yy_list)
#   })
#   agg_coords <- expand.grid(x=seq(0,x_max,length.out=nx),y=seq(0,y_max,length.out=ny))
#   coords <- agg_coords / max(x_max,y_max)
#
#   x_list <- lapply(1:num_images,\(i) {
#     matrix(1,nrow = n,ncol = p)
#   })
#
#   out <- readRDS(save_file_ltb)
#
#   # lambda <- get_rvars(out,"lambda",thin=n_thin)
#   # lambda
#   # mcmc_trace(as_draws_df(lambda))
#   # theta <- get_rvars(out,"theta",thin=n_thin)
#   # theta
#   # mcmc_trace(as_draws_df(theta[1,]))
#   #
#   hs <- seq(0,1,0.1)
#   xl <- cross_list(out,hs,thin=n_thin)
#   out_actual <- list(theta_mcmc=matrix(c(philist,rep(0,k)),nrow = 2,ncol=k,byrow=TRUE),
#                      lambda_mcmc=Lambda)
#   out_actual <- list(lapply(out_actual,\(o) {
#     dim(o) <- c(dim(o),1)
#     o
#   }))
#   #
#   # h_ix <- 5
#   # trace_df <- as_draws_df(xl[[h_ix]]) %>%
#   #   pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
#   #   separate(variable,into = c("type1","type2"),sep=",") %>%
#   #   mutate(type1 = sub("^x\\[","",type1),
#   #          type2 = sub("\\]","",type2))
#   #
#   # trace_df %>%
#   #   mutate(.chain = factor(.chain)) %>%
#   #   ggplot(aes(.iteration,value,color=.chain,group=.chain)) +
#   #   geom_line() +
#   #   facet_grid(type1~type2,
#   #              labeller = label_wrap_gen(width=8)) +
#   #   theme_bw() +
#   #   theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
#   #   theme(axis.title = element_text(size=20),
#   #         axis.text = element_text(size=12),
#   #         strip.text = element_text(size=10)) +
#   #   labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]))
#   #
#   xl_actual <- cross_list(out_actual,hs,thin=n_thin)
#   xl_actual <- lapply(xl_actual,\(x) E(x))
#
#   types <- 1:q
#   # # spatial cross-correlation over distance plots
#   # unique_combinations_with_self(types) %>%
#   #   # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
#   #   as.data.frame() %>%
#   #   as_tibble() %>%
#   #   magrittr::set_colnames(c("t1","t2")) %>%
#   #   group_by(t1,t2) %>%
#   #   group_modify(~{
#   #     ix1 <- which(types == .y$t1)
#   #     ix2 <- which(types == .y$t2)
#   #
#   #     mu <- unlist(lapply(xl,\(x) {
#   #       E(x[ix1,ix2])
#   #     }))
#   #
#   #     mu_actual <- unlist(lapply(xl_actual,\(x) {
#   #       x[ix1,ix2]
#   #     }))
#   #
#   #     diff <- unlist(lapply(xl,\(x) {
#   #       E(x[ix1,ix2])
#   #     }))
#   #
#   #     # int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value
#   #
#   #     sigma <- unlist(lapply(xl,\(x) {
#   #       sd(x[ix1,ix2])
#   #     }))
#   #
#   #     lb <- unlist(lapply(xl,\(x) {
#   #       quantile(x[ix1,ix2],probs = 0.025)
#   #     }))
#   #
#   #     ub <- unlist(lapply(xl,\(x) {
#   #       quantile(x[ix1,ix2],probs = 0.975)
#   #     }))
#   #     tibble(mu=mu,mu_actual,lb=lb,ub=ub,hs=hs)#,auc=int_val)
#   #   }) %>%
#   #   ungroup() -> xl_e
#   #
#   #
#   # xl_e %>%
#   #   mutate(hs = hs*x_max) %>%
#   #   pivot_longer(c(mu,mu_actual)) %>%
#   #   # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
#   #   # mutate(ub = mu+sigma,
#   #   #        lb = mu-sigma) %>%
#   #   ggplot() +
#   #   geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
#   #   geom_line(aes(hs,value,color=name)) +
#   #   geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
#   #   theme_minimal() +
#   #   facet_grid(t1~t2) +
#   #   theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
#   #   theme(axis.title = element_text(size=20),
#   #         axis.text = element_text(size=12),
#   #         strip.text = element_text(size=10)) +
#   #   labs(x="Distance (\u03bcm)",y="Cross-correlation")
#
#
#   unique_combinations_with_self(types) %>%
#     as.data.frame() %>%
#     as_tibble() %>%
#     magrittr::set_colnames(c("t1","t2")) %>%
#     group_by(t1,t2) %>%
#     group_modify(~{
#       ix1 <- which(types == .y$t1)
#       ix2 <- which(types == .y$t2)
#
#       mu <- as.vector(do.call(rbind,lapply(xl,\(x) {
#         x[ix1,ix2]
#       })))
#
#       ess_bulk <- as.vector(do.call(rbind,lapply(xl,\(x) {
#         posterior::ess_bulk(x[ix1,ix2])
#       })))
#
#       ess_tail <- as.vector(do.call(rbind,lapply(xl,\(x) {
#         posterior::ess_tail(x[ix1,ix2])
#       })))
#
#       rhats <- as.vector(do.call(rbind,lapply(xl,\(x) {
#         posterior::rhat(x[ix1,ix2])
#       })))
#
#       mu_actual <- unlist(lapply(xl_actual,\(x) {
#         x[ix1,ix2]
#       }))
#
#       diff <- as.vector(do.call(rbind,lapply(1:length(xl),\(i) {
#         xl[[i]][ix1,ix2] - xl_actual[[i]][ix1,ix2]
#       })))
#
#       sigma <- unlist(lapply(xl,\(x) {
#         sd(x[ix1,ix2])
#       }))
#
#       lb <- unlist(lapply(xl,\(x) {
#         quantile(x[ix1,ix2],probs = 0.025)
#       }))
#
#       ub <- unlist(lapply(xl,\(x) {
#         quantile(x[ix1,ix2],probs = 0.975)
#       }))
#       tibble(mu=mu,mu_actual,lb=lb,ub=ub,hs=hs,diff=diff,
#              ess_bulk=ess_bulk,ess_tail=ess_tail,rhat=rhats)
#     }) %>%
#     ungroup() %>%
#     mutate(sim=sim_idx) -> xl_diff
#
#   xl_diff$waic <- out[[1]]$waic
#
#   return(xl_diff)
# }) %>%
#   bind_rows() -> df_diff
#
# saveRDS(df_diff,"examples/df_diff_varyingpixels.rds")

df_diff <- readRDS("examples/df_diff_varyingpixels.rds")

grid <- grid %>%
  mutate(sim_idx=1:nrow(grid))

df_diff %>%
  rename(sim_idx=sim) %>%
  left_join(grid,by="sim_idx") -> df_diff

# paper
df_diff %>%
  mutate(abs_diff=abs(diff)) %>%
  mutate(sz=factor(sz)) %>%
  mutate(hs = hs*max(x_max,y_max)) %>%
  # mutate(mad=median(abs(diff))) %>%
  group_by(hs,sz) %>%
  # summarise(mean_diff=mean(mad))
  summarise(mean_diff=rvar_median(abs_diff)) %>%
  # print(n=nrow(.))
  # ggplot(aes(xdist=mean_diff,y=trial_k)) +
  # stat_halfeye() +
  # facet_wrap(~hs) +
  ggplot(aes(x=hs,ydist=mean_diff,color=sz,fill=sz)) +
  stat_halfeye() +
  # theme_bipps() +
  # theme(text=element_text(size=28)) +
  labs(y="MAD between true and fitted",x="Distance (\u03bcm)",color="# of pixels",fill="# of pixels")

fsave(paste0(figures_folder,"mad_distance_sz.png"))

# supplement
# WAICs are not comparable across pixel sizes, since it is a different dataset.

# supplement
df_diff %>%
  mutate(sz=factor(sz)) %>%
  mutate(hs = hs*max(x_max,y_max)) %>%
  # mutate(mad=median(abs(diff))) %>%
  group_by(hs,sz) %>%
  summarise(mean_rhat=mean(rhat,na.rm = T)) %>%
  ggplot(aes(x=hs,y=mean_rhat,color=sz,fill=sz)) +
  stat_halfeye() +
  # theme_bipps() +
  # theme(text=element_text(size=28)) +
  labs(y=TeX("Average $\\hat{R}$"),x="Distance (\u03bcm)",color="# of pixels",fill="# of pixels")
fsave(paste0(figures_folder,"rhat_sz.png"))


# supplement
df_diff %>%
  mutate(sz=factor(sz)) %>%
  mutate(hs = hs*max(x_max,y_max)) %>%
  # mutate(mad=median(abs(diff))) %>%
  group_by(hs,sz) %>%
  summarise(mean_ess=mean(ess_bulk,na.rm = T)) %>%
  ggplot(aes(x=hs,y=mean_ess,color=sz,fill=sz)) +
  stat_halfeye() +
  # theme_bipps() +
  # theme(text=element_text(size=28)) +
  labs(y="Average Bulk ESS",x="Distance (\u03bcm)",color="# of pixels",fill="# of pixels")
fsave(paste0(figures_folder,"bulk_ess_sz.png"))

# supplement
df_diff %>%
  mutate(sz=factor(sz)) %>%
  mutate(hs = hs*max(x_max,y_max)) %>%
  # mutate(mad=median(abs(diff))) %>%
  group_by(hs,sz) %>%
  summarise(mean_ess=mean(ess_tail,na.rm = T)) %>%
  ggplot(aes(x=hs,y=mean_ess,color=sz,fill=sz)) +
  stat_halfeye() +
  # theme_bipps() +
  # theme(text=element_text(size=28)) +
  labs(y="Average Tail ESS",x="Distance (\u03bcm)",color="# of pixels",fill="# of pixels")
fsave(paste0(figures_folder,"tail_ess_sz.png"))

# example convergence plots
get_sim_and_actual <- \(sim_idx) {
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
  y_list <- lapply(1:num_images,\(i) {
    yy_list <- lapply(1:q,\(j) {
      mat <- matrix(WW[[i]][,j],nrow = max_ny,ncol = max_nx)
      mat <- matrix(rpois(nrow(mat)*ncol(mat),exp(mat)),nrow=nrow(mat),ncol=ncol(mat))
      r <- raster::raster(mat)
      r_sub <- raster::aggregate(r,fact = max_nx / nx, fun = sum, expand = FALSE)
      vec_sub <- c(raster::as.matrix(r_sub))
    })
    yy <- do.call(cbind,yy_list)
  })
  agg_coords <- expand.grid(x=seq(0,x_max,length.out=nx),y=seq(0,y_max,length.out=ny))
  coords <- agg_coords / max(x_max,y_max)

  x_list <- lapply(1:num_images,\(i) {
    matrix(1,nrow = n,ncol = p)
  })

  out <- readRDS(save_file_ltb)

  hs <- seq(0,1,0.1)
  xl <- cross_list(out,hs,thin=n_thin)
  out_actual <- list(theta_mcmc=matrix(c(philist,rep(0,k)),nrow = 2,ncol=k,byrow=TRUE),
                     lambda_mcmc=Lambda)
  out_actual <- list(lapply(out_actual,\(o) {
    dim(o) <- c(dim(o),1)
    o
  }))

  xl_actual <- cross_list(out_actual,hs,thin=n_thin)
  xl_actual <- lapply(xl_actual,\(x) E(x))

  types <- 1:q
  # spatial cross-correlation over distance plots
  unique_combinations_with_self(types) %>%
    # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
    as.data.frame() %>%
    as_tibble() %>%
    magrittr::set_colnames(c("t1","t2")) %>%
    group_by(t1,t2) %>%
    group_modify(~{
      ix1 <- which(types == .y$t1)
      ix2 <- which(types == .y$t2)

      mu <- unlist(lapply(xl,\(x) {
        E(x[ix1,ix2])
      }))

      mu_actual <- unlist(lapply(xl_actual,\(x) {
        x[ix1,ix2]
      }))

      int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

      sigma <- unlist(lapply(xl,\(x) {
        sd(x[ix1,ix2])
      }))

      lb <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.025)
      }))

      ub <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.975)
      }))
      tibble(mu=mu,mu_actual,lb=lb,ub=ub,hs=hs,auc=int_val)
    }) %>%
    ungroup() -> xl_e

  list(out=out,out_actual=out_actual,xl=xl,xl_actual=xl_actual,xl_e=xl_e,y_list=y_list, WW=WW,coords=coords)
}

sim_idx <- 1
res <- get_sim_and_actual(sim_idx)
out <- res$out
out_actual <- res$out_actual
xl <- res$xl
xl_actual <- res$xl_actual
xl_e <- res$xl_e
y_list <- res$y_list
WW <- res$WW

lambda <- get_rvars(out,"lambda",thin=n_thin)
lambda
mcmc_trace(as_draws_df(lambda[,1]))
theta <- get_rvars(out,"theta",thin=n_thin)
theta
mcmc_trace(as_draws_df(theta[1,]))

h_ix <- 5
trace_df <- as_draws_df(xl[[h_ix]]) %>%
  pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
  separate(variable,into = c("type1","type2"),sep=",") %>%
  mutate(type1 = sub("^x\\[","",type1),
         type2 = sub("\\]","",type2))

hs <- seq(0,1,0.1)

# supplement?
trace_df %>%
  mutate(.chain = factor(.chain)) %>%
  mutate(type1=factor(type1,levels=1:q)) %>%
  mutate(type2=factor(type2,levels=1:q)) %>%
  group_by(type1,type2) %>%
  mutate(value = value - mean(value)) %>%
  ungroup() %>%
  ggplot(aes(.iteration,value)) +
  geom_line() +
  facet_grid(type1~type2,
             labeller = label_wrap_gen(width=8)) +
  # theme_bipps() +
  # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]*max(x_max,y_max),"\u03bcm"))
fsave(paste0(figures_folder,"trace_df_sim1_sz15.png"))


sim_idx <- 3 # do this with sim 1 and 3, to contrast convergence differences
res <- get_sim_and_actual(sim_idx)
out <- res$out
out_actual <- res$out_actual
xl <- res$xl
xl_actual <- res$xl_actual
xl_e <- res$xl_e
y_list <- res$y_list
WW <- res$WW

lambda <- get_rvars(out,"lambda",thin=n_thin)
lambda
mcmc_trace(as_draws_df(lambda[,1]))
theta <- get_rvars(out,"theta",thin=n_thin)
theta
mcmc_trace(as_draws_df(theta[1,]))

h_ix <- 5
trace_df <- as_draws_df(xl[[h_ix]]) %>%
  pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
  separate(variable,into = c("type1","type2"),sep=",") %>%
  mutate(type1 = sub("^x\\[","",type1),
         type2 = sub("\\]","",type2))

hs <- seq(0,1,0.1)

# supplement?
trace_df %>%
  mutate(.chain = factor(.chain)) %>%
  mutate(type1=factor(type1,levels=1:q)) %>%
  mutate(type2=factor(type2,levels=1:q)) %>%
  group_by(type1,type2) %>%
  mutate(value = value - mean(value)) %>%
  ungroup() %>%
  ggplot(aes(.iteration,value)) +
  geom_line() +
  facet_grid(type1~type2,
             labeller = label_wrap_gen(width=8)) +
  # theme_bipps() +
  # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]*max(x_max,y_max),"\u03bcm"))
fsave(paste0(figures_folder,"trace_df_sim3_sz45.png"))

sim_idx <- 1
res <- get_sim_and_actual(sim_idx)
xl_e <- res$xl_e

# set.seed(2025)
# filter_out <- tibble(combo=paste0(sample(1:10,6,replace=TRUE)," <--> ",paste0(sample(1:10,6,replace=TRUE))))
# supplement
dd1 <- xl_e %>%
  # mutate(combo = paste0(t1," <--> ",t2)) %>%
  # right_join(filter_out) %>%
  mutate(hs = hs*x_max) %>%
  rename(Predicted=mu,True=mu_actual) %>%
  pivot_longer(c(Predicted,True)) %>%
  mutate(pixel_size=grid$sz[sim_idx])


sim_idx <- 2
res <- get_sim_and_actual(sim_idx)
xl_e <- res$xl_e

# set.seed(2025)
# filter_out <- tibble(combo=paste0(sample(1:10,6,replace=TRUE)," <--> ",paste0(sample(1:10,6,replace=TRUE))))
# supplement
dd2 <- xl_e %>%
  # mutate(combo = paste0(t1," <--> ",t2)) %>%
  # right_join(filter_out) %>%
  mutate(hs = hs*x_max) %>%
  rename(Predicted=mu,True=mu_actual) %>%
  pivot_longer(c(Predicted,True)) %>%
  mutate(pixel_size=grid$sz[sim_idx])


sim_idx <- 3
res <- get_sim_and_actual(sim_idx)
xl_e <- res$xl_e

# set.seed(2025)
# filter_out <- tibble(combo=paste0(sample(1:10,6,replace=TRUE)," <--> ",paste0(sample(1:10,6,replace=TRUE))))
# supplement
xl_e %>%
  mutate(hs = hs*x_max) %>%
  rename(Predicted=mu,True=mu_actual) %>%
  pivot_longer(c(Predicted,True)) %>%
  mutate(pixel_size=grid$sz[sim_idx]) -> dd3

bind_rows(dd1,dd2,dd3) %>%
  mutate(pixel_size=factor(pixel_size)) %>%
  rename(Curve=name) %>%
  ggplot() +
  # geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill="grey70") +
  geom_line(aes(hs,value,color=pixel_size,linetype = Curve,group=paste(pixel_size,Curve))) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_bipps() +
  facet_grid(t1~t2) +
  # facet_wrap(~combo)
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation",color="# of pixels")
fsave(paste0(figures_folder,"cross_cor_sim_comparison.png"))

# actual on different grid sizes
sim_idx <- 1
res <- get_sim_and_actual(sim_idx)

y1 <- res$y_list[[1]][,1]
c1 <- res$coords
d1 <- cbind(c1,val=y1,im=1)

sim_idx <- 2
res <- get_sim_and_actual(sim_idx)

y2 <- res$y_list[[1]][,1]
c2 <- res$coords
d2 <- cbind(c2,val=y2,im=2)



sim_idx <- 3
res <- get_sim_and_actual(sim_idx)

y3 <- res$y_list[[1]][,1]
c3 <- res$coords
d3 <- cbind(c3,val=y3,im=3)

p1 <- d1 %>%
  ggplot(aes(x,y,fill=val)) +
  geom_tile() +
  labs(x="X",y="Y",fill="Count") +
  theme(axis.text = element_blank()) +
  theme(legend.position="top") +
  ggplot2::scale_fill_viridis_c(option="magma",    breaks = seq(floor(min(d1$val)), ceiling(max(d1$val)), by = 2),
                                labels = function(x) ifelse(x %% 2 == 0, x, "")
  )
p1

p2 <- d2 %>%
  ggplot(aes(x,y,fill=val)) +
  geom_tile() +
  labs(x="X",y="Y",fill="Count") +
  ggplot2::scale_fill_viridis_c(option="magma") +
  theme(axis.text = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position="top")


p3 <- d3 %>%
  ggplot(aes(x,y,fill=val)) +
  geom_tile() +
  labs(x="X",y="Y",fill="Count") +
  theme(axis.text = element_blank(),
        axis.title.y = element_blank()) +
  theme(legend.position="top") +
  ggplot2::scale_fill_viridis_c(option="magma")


p1 + p2 + p3 + plot_annotation(tag_levels = 'a')
fsave(paste0(figures_folder,"Y_list_diff_res.png"))



# difference between groups (within same combination), low-res
idx1 <- 1
idx2 <- 4

res1 <- get_sim_and_actual(idx1)
out1 <- res1$out
out1_actual <- res1$out_actual
xl1 <- res1$xl
xl1_actual <- res1$xl_actual
xl1_e <- res1$xl_e

res2 <- get_sim_and_actual(idx2)
out2 <- res2$out
out2_actual <- res2$out_actual
xl2 <- res2$xl
xl2_actual <- res2$xl_actual
xl2_e <- res2$xl_e

# supplement
xl1_e %>%
  mutate(hs = hs*x_max) %>%
  pivot_longer(c(mu,mu_actual)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_minimal() +
  facet_grid(t1~t2) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  ggtitle("Group 1")

xl2_e %>%
  mutate(hs = hs*x_max) %>%
  pivot_longer(c(mu,mu_actual)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_minimal() +
  facet_grid(t1~t2) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  ggtitle("Group 2")

group_diff_actual <- lapply(1:length(xl1_actual) ,\(i){
  xl1_actual[[i]] - xl2_actual[[i]]
})

group_diff_est <- lapply(1:length(xl1) ,\(i){
  xl1[[i]] - xl2[[i]]
})

types <- 1:q

unique_combinations_with_self(types) %>%
  # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types == .y$t1)
    ix2 <- which(types == .y$t2)

    mu <- unlist(lapply(group_diff_est,\(x) {
      E(x[ix1,ix2])
    }))

    mu_actual <- unlist(lapply(group_diff_actual,\(x) {
      x[ix1,ix2]
    }))


    lb <- unlist(lapply(group_diff_est,\(x) {
      quantile(x[ix1,ix2],probs = 0.025)
    }))

    ub <- unlist(lapply(group_diff_est,\(x) {
      quantile(x[ix1,ix2],probs = 0.975)
    }))
    tibble(diff=mu,diff_actual=mu_actual,lb=lb,ub=ub,hs=hs)
  }) %>%
  ungroup() -> group_diff_e

# paper
group_diff_e %>%
  filter(t1 %in% c(1),t2%in% 1:3) %>%
  mutate(hs = hs*x_max) %>%
  rename(Predicted=diff,Observed=diff_actual) %>%
  pivot_longer(c(Predicted,Observed)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_bipps() +
  # facet_grid(t1~t2) +
  facet_wrap(~t1 + t2,ncol=5) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation",color="Difference")
  # theme(strip.background = element_blank(),
        # strip.text.x = element_blank())
# ggtitle("Difference between group 1 and group 2")
fsave(paste0(figures_folder,"group_diff_sz_lo-res.png"))

# difference between groups (within same combination), hi-res
idx1 <- 3
idx2 <- 6

res1 <- get_sim_and_actual(idx1)
out1 <- res1$out
out1_actual <- res1$out_actual
xl1 <- res1$xl
xl1_actual <- res1$xl_actual
xl1_e <- res1$xl_e

res2 <- get_sim_and_actual(idx2)
out2 <- res2$out
out2_actual <- res2$out_actual
xl2 <- res2$xl
xl2_actual <- res2$xl_actual
xl2_e <- res2$xl_e

# supplement
xl1_e %>%
  mutate(hs = hs*x_max) %>%
  pivot_longer(c(mu,mu_actual)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_minimal() +
  facet_grid(t1~t2) +
  # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  ggtitle("Group 1")

xl2_e %>%
  mutate(hs = hs*x_max) %>%
  pivot_longer(c(mu,mu_actual)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_minimal() +
  facet_grid(t1~t2) +
  # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  ggtitle("Group 2")

group_diff_actual <- lapply(1:length(xl1_actual) ,\(i){
  xl1_actual[[i]] - xl2_actual[[i]]
})

group_diff_est <- lapply(1:length(xl1) ,\(i){
  xl1[[i]] - xl2[[i]]
})

types <- 1:q

unique_combinations_with_self(types) %>%
  # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types == .y$t1)
    ix2 <- which(types == .y$t2)

    mu <- unlist(lapply(group_diff_est,\(x) {
      E(x[ix1,ix2])
    }))

    mu_actual <- unlist(lapply(group_diff_actual,\(x) {
      x[ix1,ix2]
    }))


    lb <- unlist(lapply(group_diff_est,\(x) {
      quantile(x[ix1,ix2],probs = 0.025)
    }))

    ub <- unlist(lapply(group_diff_est,\(x) {
      quantile(x[ix1,ix2],probs = 0.975)
    }))
    tibble(diff=mu,diff_actual=mu_actual,lb=lb,ub=ub,hs=hs)
  }) %>%
  ungroup() -> group_diff_e

# paper
group_diff_e %>%
  filter(t1 %in% c(1),t2%in% 1:3) %>%
  mutate(hs = hs*x_max) %>%
  rename(Predicted=diff,Observed=diff_actual) %>%
  pivot_longer(c(Predicted,Observed)) %>%
  ggplot() +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line(aes(hs,value,color=name)) +
  geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
  # theme_bipps() +
  # facet_grid(t1~t2) +
  facet_wrap(~t1+t2) +
  # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
        # axis.text = element_text(size=12),
        # strip.text = element_text(size=10),
        # strip.background = element_blank(),
        # strip.text.x = element_blank()) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation",color="Difference")
# ggtitle("Difference between group 1 and group 2")
fsave(paste0(figures_folder,"group_diff_sz_hi-res.png"))
