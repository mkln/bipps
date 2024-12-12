devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(bayesplot)

args <- commandArgs(trailingOnly = TRUE)
sim_idx <- as.integer(args[1])

# simset1
actual_ks <- 3
trial_ks <- c(2,3,6)

# simset2
# actual_ks <- 5
# trial_ks <- c(3,5,7)

# mcmc and model settings
grid <- expand.grid(actual_k=actual_ks,trial_k=trial_ks,sim=1:10)
start_idx <- 1
file_prefix <- "/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset1_group_diff"
seed_start <- 24

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
mu <- -2
q <- 10

# simulation settings
nx <- 30
ny <- 30
n <- nx*ny
x_max <- 1919
y_max <- 1439
sigmasq <- 1

p <- 1

coords <- expand.grid(x=seq(0,x_max,length.out=nx),y=seq(0,y_max,length.out=ny))

coords <- coords / max(x_max,y_max)
c_mat <- as.matrix(coords)
d_coords <- as.matrix(dist(c_mat))

print(sim_idx)
seed <- seed_start + sim_idx
actual_k <- grid$actual_k[sim_idx]
trial_k <- grid$trial_k[sim_idx]

set.seed(seed)

save_file <- paste0(file_prefix,"_",sim_idx,".rds")
save_file_ltb <- paste0(file_prefix,"_",sim_idx,"_ltb.rds")

philist <- runif(actual_k,phi_range[1],phi_range[2])

LClist <- 1:actual_k %>% lapply(\(i) t(chol(
  exp(- philist[i] * d_coords)
)))

VV <- lapply(1:num_images,\(j) {
  wlist <- lapply(1:actual_k,\(i) LClist[[i]] %*% rnorm(n))

  # factor matrix
  do.call(cbind, wlist)
})


# factor loadings
Lambda <- matrix(0, q, actual_k)
diag(Lambda) <- runif(actual_k, 0.5, 1)
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

y_list <- lapply(WW,\(ww) {
  matrix(rpois(nrow(ww)*ncol(ww),exp(ww)),nrow=nrow(ww),ncol=ncol(ww))
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
                                        k = trial_k,
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

out_ltb <- lapply(out,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc,beta_mcmc=o$beta_mcmc))

saveRDS(out_ltb,save_file_ltb)

if(do_plots) {
  out <- readRDS("out_simchol_lt.rds")

  lambda <- get_rvars(out,"lambda",thin=n_thin)
  lambda
  mcmc_trace(as_draws_df(lambda))
  theta <- get_rvars(out,"theta",thin=n_thin)
  theta
  mcmc_trace(as_draws_df(theta[1,]))

  hs <- seq(0,1,0.1)
  xl <- cross_list(out,hs,thin=n_thin)
  out_actual <- list(theta_mcmc=matrix(c(philist,rep(0,k)),nrow = 2,ncol=k,byrow=TRUE),
                     lambda_mcmc=Lambda)
  out_actual <- list(lapply(out_actual,\(o) {
    dim(o) <- c(dim(o),1)
    o
  }))

  h_ix <- 5
  trace_df <- as_draws_df(xl[[h_ix]]) %>%
    pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
    separate(variable,into = c("type1","type2"),sep=",") %>%
    mutate(type1 = sub("^x\\[","",type1),
           type2 = sub("\\]","",type2))

  trace_df %>%
    mutate(.chain = factor(.chain)) %>%
    ggplot(aes(.iteration,value,color=.chain,group=.chain)) +
    geom_line() +
    facet_grid(type1~type2,
               labeller = label_wrap_gen(width=8)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=12),
          strip.text = element_text(size=10)) +
    labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]))

  xl_actual <- cross_list(out_actual,hs,thin=n_thin)
  xl_actual <- lapply(xl_actual,\(x) E(x))

  unique_combinations_with_self <- function(data) {
    # Generate all pairs including self-pairings
    all_pairs <- expand.grid(data, data)

    # Sort the pairs and remove duplicates
    all_pairs <- t(apply(all_pairs, 1, sort))
    unique_pairs <- unique(all_pairs)

    return(unique_pairs)
  }

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

      # int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

      sigma <- unlist(lapply(xl,\(x) {
        sd(x[ix1,ix2])
      }))

      lb <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.025)
      }))

      ub <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.975)
      }))
      tibble(mu=mu,mu_actual,lb=lb,ub=ub,hs=hs)#,auc=int_val)
    }) %>%
    ungroup() -> xl_e


  xl_e %>%
    mutate(hs = hs*x_max) %>%
    pivot_longer(c(mu,mu_actual)) %>%
    # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
    # mutate(ub = mu+sigma,
    #        lb = mu-sigma) %>%
    ggplot() +
    geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
    geom_line(aes(hs,value,color=name)) +
    geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
    theme_minimal() +
    facet_grid(t1~t2) +
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=12),
          strip.text = element_text(size=10)) +
    labs(x="Distance (\u03bcm)",y="Cross-correlation")
}
