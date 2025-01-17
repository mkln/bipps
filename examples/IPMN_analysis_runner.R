devtools::load_all()
library(spatstat)
library(tidyverse)

set.seed(2020)

n_samples <- 1000
n_burnin <- 20000
n_thin <- 10
n_threads <- 8
block_size <- 50
ks <- c(2,4,6)
starting <- list(phi = 5)
prior <- list(phi = c(0.1,10))
save_file <- "IPMN_varyingk_nburn20k_nthin10_nsamp1e3_chain1.rds"
save_file_lt <- "IPMN_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds"

chains <- 1

dat1 <- readRDS("examples/data/PDAC_POS_DATA.rds") %>%
  rename(X=Cell.X.Position,
         Y=Cell.Y.Position,
         type=Cellname,
         Spot=SlideID)

dat2 <- readRDS("examples/data/IPMN_POS_DATA.rds") %>%
  rename(X=Cell.X.Position,
         Y=Cell.Y.Position,
         type=Cellname,
         Spot=SlideID)

# dat1 %>%
#   filter(Spot == "PATID_538_SLIDE_1") %>%
#   ggplot(aes(X,Y,color=type)) +
#   geom_point() +
#   scale_color_manual(values=as.vector(pals::glasbey())) +
#   theme_bw()

dat1 %>%
  distinct(Spot) %>%
  filter(!(Spot %in% c("PATID_12_SLIDE_1","PATID_538_SLIDE_1"))) %>% # these are bad images
  pull(Spot) -> spots1

dat2 %>%
  distinct(Spot) %>%
  pull(Spot) -> spots2

dat1 <- dat1 %>%
  filter(Spot %in% spots1)

dat2 <- dat2 %>%
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

out1 <- pixellate_grid(dat1$X,dat1$Y,dat1$type,dat1$Spot,nx,ny)
y_list1 <- out1$y_list
coords1 <- out1$coords
x_list1 <- lapply(y_list1,\(yy) {
  matrix(1,nrow = nrow(yy),ncol = 1)
})

out2 <- pixellate_grid(dat2$X,dat2$Y,dat2$type,dat2$Spot,nx,ny)
y_list2 <- out2$y_list
coords2 <- out2$coords
x_list2 <- lapply(y_list2,\(yy) {
  matrix(1,nrow = nrow(yy),ncol = 1)
})

# p <- plot_y_list(y_list2,coords2)
# p[[52]]

out2 <- lapply(ks,\(k) multi_bipps(y_list2,
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
                                   just_preprocess = F))
saveRDS(out2,save_file)
out2_lt <- lapply(out2,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc,waic=o$waic,beta_mcmc=o$beta_mcmc))
saveRDS(out2_lt,save_file_lt)
