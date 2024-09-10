rm(list = ls())
devtools::load_all()
library(spatstat)
library(tidyverse)

set.seed(2020)

df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  # dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)
# mutate(Spot = factor(Spot))


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

nx <- ny <- 20
out1 <- create_y_list(dat1$X,dat1$Y,dat1$type,dat1$Spot,nx,ny)
y_list1 <- out1$y_list
coords1 <- out1$coords
x_list1 <- lapply(y_list1,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

# p <- plot_y_list(y_list1,coords1)
# p[[49]]

out2 <- create_y_list(dat2$X,dat2$Y,dat2$type,dat2$Spot,nx,ny)
y_list2 <- out2$y_list
coords2 <- out2$coords
x_list2 <- lapply(y_list2,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

n_samples <- 2000
n_burnin <- 1000
n_thin <- 1
n_threads <- 16
block_size <- 35
k <- 4
starting <- list(phi = 100)
prior <- list(phi = c(0.1, 200))

out1 <- multi_bipps(y_list1,
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
                      verbose = F, debug = F
                    ),
                    just_preprocess = F)
saveRDS(out1,"out1_CRC_analysis.rds")

out2 <- multi_bipps(y_list2,
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
                      verbose = F, debug = F
                    ),
                    just_preprocess = F)
saveRDS(out2,"out2_CRC_analysis.rds")

out1 <- readRDS("out1_CRC_analysis.rds")
out2 <- readRDS("out2_CRC_analysis.rds")

hs <- seq(0,0.5,0.05)

xl1 <- cross_list(out1,hs)
xl2 <- cross_list(out2,hs)

xl_diff <- lapply(1:length(hs),\(i) {
  x1 <- xl1[[i]]
  x2 <- xl2[[i]]

  x1 - x2
})

unique_combinations_with_self <- function(data) {
  # Generate all pairs including self-pairings
  all_pairs <- expand.grid(data, data)

  # Sort the pairs and remove duplicates
  all_pairs <- t(apply(all_pairs, 1, sort))
  unique_pairs <- unique(all_pairs)

  return(unique_pairs)
}

library(posterior)
unique_combinations_with_self(types_intersect) %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types_intersect == .y$t1)
    ix2 <- which(types_intersect == .y$t2)

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
