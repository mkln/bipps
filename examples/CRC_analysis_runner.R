rm(list = ls())
devtools::load_all()
library(spatstat)
library(tidyverse)

set.seed(2020)

df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  # dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)
# mutate(Spot = factor(Spot))

# df_raw %>%
#   filter(Spot == "42_A") %>%
#   ggplot(aes(X,Y,color=type)) +
#   geom_point() +
#   scale_color_manual(values=as.vector(pals::glasbey())) +
#   theme_bw()

df_raw %>%
  distinct(Spot,groups) %>%
  filter(groups == 1) %>%
  filter(!(Spot %in% c("67_B","57_A"))) %>% # these are really weird images
  pull(Spot) -> spots1

df_raw %>%
  distinct(Spot,groups) %>%
  filter(groups == 2) %>%
  pull(Spot) -> spots2

dat1 <- df_raw %>%
  filter(Spot %in% spots1)

dat2 <- df_raw %>%
  filter(Spot %in% spots2)

dat1 %>%
  count(Spot,type) %>%
  group_by(type) %>%
  summarise(mn = median(n)) %>%
  filter(mn > 10) %>%
  pull(type) -> types1

dat2 %>%
  count(Spot,type) %>%
  group_by(type) %>%
  summarise(mn = median(n)) %>%
  filter(mn > 10) %>%
  pull(type) -> types2


types_intersect <- intersect(types1,types2)

dat1 <- dat1 %>%
  filter(type %in% types_intersect)

dat2 <- dat2 %>%
  filter(type %in% types_intersect)

bind_rows(dat1,dat2) %>%
  dplyr::group_by(Spot) %>%
  dplyr::mutate(X = X - min(X),
                Y = Y - min(Y)) %>%
  dplyr::ungroup() %>%
  select(X,Y) %>%
  apply(.,2,max) -> max_dim

# x <- dat1$X
# y <- dat1$Y
# types <- dat1$type
# image_ids <- dat1$Spot

pix_dim <- 70
nx <- ceiling(max_dim[1]/pix_dim)
ny <- ceiling(max_dim[2]/pix_dim)
# nx <- 20
# ny <- 20

out1 <- pixellate_grid(dat1$X,dat1$Y,dat1$type,dat1$Spot,nx,ny)
y_list1 <- out1$y_list
coords1 <- out1$coords
x_list1 <- lapply(y_list1,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

# p <- plot_y_list(y_list1,coords1)
# p[[47]]

out2 <- pixellate_grid(dat2$X,dat2$Y,dat2$type,dat2$Spot,nx,ny)
y_list2 <- out2$y_list
coords2 <- out2$coords
x_list2 <- lapply(y_list2,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

n_samples <- 1000
n_burnin <- 10000
n_thin <- 1
n_threads <- 16
block_size <- 50
k <- 4
starting <- list(phi = 5)
# D <- sqrt(max_dim[1]^2 + max_dim[2]^2)
# a <- -log(0.05) * 3 / D
prior <- list(phi = c(0.1,10))

chains <- 4

out1 <- lapply(1:chains,\(i) multi_bipps(y_list1,
                    x_list1,
                    coords1,
                    k = k,
                    family = "poisson",
                    block_size = block_size,
                    n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                    n_threads = n_threads,
                    starting = starting,
                    prior = prior,
                    settings = list(adapting = T, saving = T, ps = T),
                    verbose = 10,
                    debug = list(
                      sample_beta = T, sample_tausq = F,
                      sample_theta = T, sample_w = T, sample_lambda = T,
                      verbose = T, debug = F
                    ),
                    just_preprocess = F))
saveRDS(out1,"out1_chains_CRC_analysis.rds")

out2 <- lapply(1:chains,\(i) multi_bipps(y_list2,
                    x_list2,
                    coords2,
                    k = k,
                    family = "poisson",
                    block_size = block_size,
                    n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                    n_threads = n_threads,
                    starting = starting,
                    prior = prior,
                    settings = list(adapting = T, saving = T, ps = T),
                    verbose = 10,
                    debug = list(
                      sample_beta = T, sample_tausq = F,
                      sample_theta = T, sample_w = T, sample_lambda = T,
                      verbose = T, debug = F
                    ),
                    just_preprocess = F))
saveRDS(out2,"out2_chains_CRC_analysis.rds")
