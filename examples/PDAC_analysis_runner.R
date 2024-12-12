devtools::load_all()
library(spatstat)
library(tidyverse)

set.seed(2020)

n_samples <- 1000
n_burnin <- 5000
n_thin <- 1
n_threads <- 16
block_size <- 50
k <- 2
starting <- list(phi = 5)
# D <- sqrt(max_dim[1]^2 + max_dim[2]^2)
# a <- -log(0.05) * 3 / D
prior <- list(phi = c(0.1,10))
save_file <- "PDAC_analysis.rds"
save_file_lt <- "PDAC_analysis_lt.rds"

chains <- 2

pdac_raw <- readRDS("examples/data/PDAC_DATA.rds")
  # dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)
