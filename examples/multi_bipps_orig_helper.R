rm(list = ls())
devtools::load_all()
library(magrittr)
library(spatstat)
library(pals)
library(tidyverse)

set.seed(2020)

# df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
#   dplyr::mutate(type = as.factor(type)) %>%
#   dplyr::rename(Spot = spots)
#
# dat <- df_raw %>%
#   filter(Spot %in% c("59_A")) %>%
#   droplevels()
  #,"59_B","60_A","60_B"))
#
# rm(df_raw)

# write_csv(dat,"examples/data/dat_59_A.csv")
dat <- read_csv("examples/data/dat_59_A.csv") %>%
  mutate(type = as.factor(type))

# all.equal(dat,dat2)


x <- dat$X
y <- dat$Y
types <- dat$type
image_ids <- dat$Spot
nx <- ny <- 20
covariates <- NULL
family <- "poisson"
n_partition <- 8
k <- NULL
verbose <- 20
n_samples <- 100
n_burnin <- 100
n_thin <- 1

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
             verbose=FALSE, debug=FALSE)
n_threads <- 1

out <- create_y_list(x,y,types,image_ids,nx,ny)
y_list <- out$y_list
coords <- out$grid %>%
  ungroup() %>%
  select(x,y) %>%
  as.matrix()

x_list <- lapply(1:length(unique(image_ids)),\(i) {
  matrix(0,nrow=nrow(coords),ncol=1)
})

# multi_bipps(y_list,x_list,coords,verbose=10,prior=prior,n_threads=1)
ncol(y_list[[1]])

# y_list <- lapply(y_list,\(yy) {
#   yy[is.na(yy)] <- 0
#   yy
# })
#
# y_list <- lapply(y_list,\(yy) {
#   colnames(yy) <- NULL
#   yy
# })

# list_args <- list(
#   y_list=y_list, x_list=x_list, coords=coords,
#   k = 5,
#   family = "poisson",
#   block_size = 5,
#   n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
#   n_threads = 1,
#   starting = list(phi = 1),
#   prior = list(phi = c(1, 200), nu = c(.5, .5)),
#   settings = list(adapting = T, saving = F, ps = T, hmc = 0),
#   verbose = 10,
#   debug = list(
#     sample_beta = F, sample_tausq = F,
#     sample_theta = T, sample_w = T, sample_lambda = T,
#     verbose = T, debug = T
#   )
# )

# saveRDS(list_args,"list_args_df_raw.rds")
#
# list_args_dat <- readRDS("list_args_dat.rds")
# list_args_df_raw <- readRDS("list_args_df_raw.rds")
#
# # list_args_dat$y_list[[1]][3,3] <- NA
#
# all.equal(list_args_dat,list_args_df_raw)
#
# # ok so these have exactly the same arguments coming in. Are there some global variables causing havoc?
# library(codetools)
# findGlobals(multi_bipps, merge=FALSE)[['variables']]

# there are some globals here. may want to fix those.

meshout <- multi_bipps(y_list, x_list, coords,
                 k = 5,
                 family = "poisson",
                 block_size = 5,
                 n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                 n_threads = 1,
                 starting = list(phi = 1),
                 prior = list(phi = c(1, 200), nu = c(.5, .5)),
                 settings = list(adapting = T, saving = T, ps = T, hmc = 0),
                 verbose = 10,
                 debug = list(
                   sample_beta = F, sample_tausq = F,
                   sample_theta = T, sample_w = T, sample_lambda = T,
                   verbose = T, debug = T
                 ),
                 just_preprocess = T
)

saveRDS(meshout,"meshout_main.rds")
# meshout_dat <- readRDS("meshout_dat.rds")
# meshout_df_raw <- readRDS("meshout_df_raw.rds")
#
# all.equal(meshout_dat,meshout_df_raw)
