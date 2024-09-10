rm(list = ls())
devtools::load_all()
library(magrittr)
library(spatstat)
library(pals)
library(tidyverse)
library(tidybayes)

set.seed(2020)

df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  # dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)
  # mutate(Spot = factor(Spot))

# 59_A,B are in group 2
#
df_raw %>%
  distinct(Spot,groups)
#
dat2 <- df_raw %>%
  filter(Spot %in% c("59_A","59_B"))

dat1 <- df_raw %>%
  filter(Spot %in% c("1_A","1_B"))

types1 <- unique(dat1$type)
types2 <- unique(dat2$type)

types_intersect <- intersect(types1,types2)

dat1 <- dat1 %>%
  filter(type %in% types_intersect)

dat2 <- dat2 %>%
  filter(type %in% types_intersect)
  # group_by(Spot) %>%
  # mutate(type = fct_lump_n(type,n=20)) %>%
  # filter(type != "Other") %>%
  # droplevels()

# dat <- dat %>%
#   group_by(type) %>%
#   filter(n_distinct(Spot) == n_distinct(dat$Spot)) %>%
#   ungroup() %>%
#   droplevels()
# dat <- read_csv("examples/data/dat_59_A.csv")
#   mutate(type = as.factor(type))

nx <- ny <- 20
verbose <- 20
out1 <- create_y_list(dat1$X,dat1$Y,dat1$type,dat1$Spot,nx,ny)
y_list1 <- out1$y_list
coords1 <- out1$coords
# p <- plot_y_list(y_list[1],coords)
# p
x_list1 <- lapply(y_list1,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

out2 <- create_y_list(dat2$X,dat2$Y,dat2$type,dat2$Spot,nx,ny)
y_list2 <- out2$y_list
coords2 <- out2$coords
# p <- plot_y_list(y_list[1],coords)
# p
x_list2 <- lapply(y_list2,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

# covariates <- NULL
family <- "poisson"
# n_partition <- 8
k <- 4
n_samples <- 2000
n_burnin <- 1000
n_thin <- 1
#
settings = list(adapting=TRUE,
                ps=TRUE, saving=TRUE, low_mem=FALSE, hmc=0)
prior = list(beta=NULL, tausq=NULL, sigmasq = NULL,
             phi=c(1, 40), a=NULL, nu = NULL,
             toplim = NULL, btmlim = NULL, set_unif_bounds=NULL)
starting = list(beta=NULL, tausq=NULL, theta=NULL,
                lambda=NULL, v=NULL,  a=NULL, nu = NULL,
                mcmcsd=.05, mcmc_startfrom=0)
debug = list(sample_beta=TRUE, sample_tausq=TRUE,
             sample_theta=TRUE, sample_w=TRUE, sample_lambda=TRUE,
             verbose=TRUE, debug=TRUE)
n_threads <- 3
axis_partition <- NULL
block_size <- 25


out1 <- multi_bipps(y_list1,
                   x_list1,
                   coords1,
                   k = 4,
                   family = "poisson",
                   block_size = 25,
                   n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                   n_threads = n_threads,
                   starting = list(phi = 100),
                   prior = list(phi = c(0.1, 200)),
                   settings = list(adapting = T, saving = T, ps = T),
                   verbose = 10,
                   debug = list(
                     sample_beta = T, sample_tausq = F,
                     sample_theta = T, sample_w = T, sample_lambda = T,
                     verbose = F, debug = F
                   ),
                   just_preprocess = F)

out2 <- multi_bipps(y_list2,
                    x_list2,
                    coords2,
                    k = 4,
                    family = "poisson",
                    block_size = 25,
                    n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                    n_threads = n_threads,
                    starting = list(phi = 100),
                    prior = list(phi = c(0.1, 200)),
                    settings = list(adapting = T, saving = T, ps = T),
                    verbose = 10,
                    debug = list(
                      sample_beta = T, sample_tausq = F,
                      sample_theta = T, sample_w = T, sample_lambda = T,
                      verbose = F, debug = F
                    ),
                    just_preprocess = F)
# lambda <- get_rvars(out,"lambda")
# corrs1 <- get_cross_cor(out1,0.1)
# theta1 <- get_rvars(out,"theta")

# also look at correlations over space as well - use first row of theta
# find distance for which correlation becomes low enough

# solve following function, for each draw (identify where equal to 0.05)
# function(h){
#
# }

hs <- seq(0,0.5,0.05)

xl1 <- cross_list(out1,hs)
xl2 <- cross_list(out2,hs)

xl_diff <- lapply(1:length(hs),\(i) {
  x1 <- xl1[[i]]
  x2 <- xl2[[i]]

  x1 - x2
})

types <- colnames(y_list[[1]])
unique_combinations_with_self <- function(data) {
  # Generate all pairs including self-pairings
  all_pairs <- expand.grid(data, data)

  # Sort the pairs and remove duplicates
  all_pairs <- t(apply(all_pairs, 1, sort))
  unique_pairs <- unique(all_pairs)

  return(unique_pairs)
}

unique_combinations_with_self(types) %>%
  as.data.frame() %>%
  as_tibble() %>%
  set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types == .y$t1)
    ix2 <- which(types == .y$t2)

    mu <- unlist(lapply(xl_diff,\(x) {
      E(x[ix1,ix2])
    }))

    int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

    sigma <- unlist(lapply(xl_diff,\(x) {
      sd(x[ix1,ix2])
    }))

    lb <- unlist(lapply(xl_diff,\(x) {
      quantile(x[ix1,ix2],probs = 0.025)
    }))

    ub <- unlist(lapply(xl_diff,\(x) {
      quantile(x[ix1,ix2],probs = 0.975)
    }))
    tibble(mu=mu,lb=lb,ub=ub,hs=hs,auc=int_val)
  }) %>%
  ungroup() -> xldiff_e

xldiff_e %>%
  filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red") +
  theme_bw() +
  facet_wrap(t1~t2) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))

rownames(corrs) <- colnames(corrs) <- colnames(y_list[[1]])

co

upper_corrs <- corrs
upper_corrs[lower.tri(upper_corrs)] <- NA

summarise_draws(upper_corrs) %>%
  separate(variable,sep=",",into=c("x","y")) %>%
  mutate(y=factor(gsub("]","",y),levels=rownames(corrs)),
         x=factor(gsub("upper_corrs\\[","",x),levels=rownames(corrs))) %>%
  # filter(x != y) %>%
  filter(!is.na(mean)) %>%
  mutate(mean = ifelse(x==y,NA,mean)) %>%
  ggplot(aes(x,y,fill=mean)) +
  geom_tile() +
  scico::scale_fill_scico(palette = "cork",direction = 1,midpoint=0) +
  labs(x="Type 1",y="Type 2") +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))


library(bayesplot)
library(posterior)

# colnames(lambda) <- rownames(lambda) <- colnames(y_list[[1]])
summarise_draws(lambda) %>%
  mutate(variable = gsub("lambda","",variable)) %>%
  print(n=nrow(.))

summarise_draws(theta[1,])
mcmc_areas_data(as_draws_df(lambda))
mcmc_dens(as_draws_df(lambda))
mcmc_trace(as_draws_df(theta[1,]))

mcmc_trace(as_draws_df(corrs[1:3,1:3]))

mcmc_pairs(as_draws_df(corrs[1:3,1:3]))



phi <- theta[1,] %>%
  as_draws_matrix()

colnames(phi) <-paste0("X_",1:15)

mcmc_pairs(as_draws_rvars(theta[1,1:5]))

# saveRDS(out,"meshout_cleanup_orig.rds")
# #
# out_cleanup <- readRDS("meshout_cleanup_orig.rds")
# out_main <- readRDS("meshout_main.rds")
# all.equal(out_cleanup,out_main)
## assuming that domain for all images has the same shape (scaling factor) - could be maximum over all input images.
# shift to bottom-left corner (where DAG begins).
# shrinking or expanding images - need an offset probably.
# vs adding pixels of the same size to smaller images.
# for now - use biggest image as domain, put the smaller images in the bottom left of the domain. boundary areas should be NA or zero.
# alt: find nearest grid cell in pre-constructed grid.
#

