devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(tidybayes)
library(bayesplot)
library(patchwork)

theme_set(theme_bw(base_size=11, base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))

# latex_width <- 6
fsave <- \(fname) {
  # if(ar == "flat") {
  #   height = latex_width * 0.66
  # } else if(ar == "square") {
  #   height = latex_width
  # }
  ggsave(paste0(figures_folder,fname),dpi=300, height=5, width=8, units="in")
}
figures_folder <- "examples/data/figures/CRC_postpred/"
ks <- c(2,4,6,8,10)

n_thin <- 10
df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  # dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)
# mutate(Spot = factor(Spot))


# pp plot
df <- df_raw %>%
  filter(Spot == "59_A") %>%
  mutate(type = fct_lump_min(type,min=20)) %>%
  filter(type != "Other") %>%
  droplevels() %>%
  rename(Type=type)

df_raw %>%
  dplyr::group_by(Spot) %>%
  dplyr::mutate(X = X - min(X),
                Y = Y - min(Y)) %>%
  dplyr::ungroup() %>%
  select(X,Y) %>%
  apply(.,2,max) -> max_dim


pix_dim <- 70
nx <- ceiling(max_dim[1]/pix_dim)
ny <- ceiling(max_dim[2]/pix_dim)

out <- pixellate_grid(df$X,df$Y,df$Type,df$Spot,nx,ny)

y_list <- out$y_list
coords <- out$coords

p1 <- plot_y_list(y_list,coords)
p1
# 54_B
# 53_B

max_range <- max(max_dim)

# model fitting
# group 1 is CLR
# group 2 is DII
df_raw %>%
  distinct(Spot,groups) %>%
  filter(groups == 1) %>%
  filter(!(Spot %in% c("67_B","57_A"))) %>% # these are really weird images
  pull(Spot) -> spots1

df_raw %>%
  distinct(Spot,groups) %>%
  filter(groups == 2) %>%
  filter(!(Spot %in% c("53_B","54_B"))) %>% # these are also really weird images
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

k_idx <- 3 # corresponds to k = 6

out1 <- readRDS("examples/data/group1_CRC_intercept_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")[k_idx]
out2 <- readRDS("examples/data/group2_CRC_intercept_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")[k_idx]

# generate potential subjects
seed <- 2025
num_images <- 1
k <- ks[k_idx]
c_mat <- as.matrix(coords)
d_coords <- as.matrix(dist(c_mat))
n <- nrow(coords)
p <- 1

phi1 <- as.vector(E(get_rvars(out1,"theta")[1,]))
phi2 <- as.vector(E(get_rvars(out2,"theta")[1,]))

LClist1 <- 1:k %>% lapply(\(i) t(chol(
  exp(- phi1[i] * d_coords)
)))

LClist2 <- 1:k %>% lapply(\(i) t(chol(
  exp(- phi2[i] * d_coords)
)))

set.seed(seed)
z_var <- lapply(1:num_images, \(i) rnorm(n))

VV1 <- lapply(1:num_images,\(j) {
  wlist <- lapply(1:k,\(i) LClist1[[i]] %*% z_var[[j]])

  # factor matrix
  do.call(cbind, wlist)
})

VV2 <- lapply(1:num_images,\(j) {
  wlist <- lapply(1:k,\(i) LClist2[[i]] %*% z_var[[j]])

  # factor matrix
  do.call(cbind, wlist)
})

Lambda1 <- E(get_rvars(out1,"lambda"))
Lambda2 <- E(get_rvars(out2,"lambda"))

Beta1 <- E(get_rvars(out1,"beta"))
Beta2 <- E(get_rvars(out2,"beta"))

x_list <- lapply(1:num_images,\(i) {
  matrix(1,nrow = n,ncol = p)
})

WW1 <- lapply(1:num_images,\(i) {
  mat <- VV1[[i]] %*% t(Lambda1) + x_list[[i]] %*% Beta1
  # print(range(exp(mat)))
  colnames(mat) <- types_intersect
  mat
})

WW2 <- lapply(1:num_images,\(i) {
  mat <- VV2[[i]] %*% t(Lambda2) + x_list[[i]] %*% Beta2
  # print(range(exp(mat)))
  colnames(mat) <- types_intersect
  mat
})

set.seed(seed)
y_list1 <- lapply(WW1,\(ww) {
  mat <- matrix(rpois(nrow(ww)*ncol(ww),exp(ww)),nrow=nrow(ww),ncol=ncol(ww))
  colnames(mat) <- types_intersect
  mat
})

set.seed(seed)
y_list2 <- lapply(WW2,\(ww) {
  mat <- matrix(rpois(nrow(ww)*ncol(ww),exp(ww)),nrow=nrow(ww),ncol=ncol(ww))
  colnames(mat) <- types_intersect
  mat
})

# p1 <- plot_y_list(y_list1,coords)
# p1[[1]]
#
# p2 <- plot_y_list(y_list2,coords)
# p2[[1]]
#
# # also plot latent factors
# p1 <- plot_y_list(WW1,coords)
# p1[[1]]
#
# p2 <- plot_y_list(WW2,coords)
# p2[[1]]

compare_postpred <- function(y_list1,y_list2,coords) {
  plots <- lapply(1:length(y_list1),\(i) {
    y1 <- y_list1[[i]]
    y2 <- y_list2[[i]]
    y1 <- dplyr::bind_cols(y1,coords)
    y2 <- dplyr::bind_cols(y2,coords)
    y1$group <- 1
    y2$group <- 2
    dplyr::bind_rows(y1,y2) %>%
      tidyr::pivot_longer(-c(x,y,group),names_to = "type",values_to = "count") %>%
      dplyr::mutate(x=x*max_range,
                    y=y*max_range) %>%
      ggplot2::ggplot(ggplot2::aes(x,y,fill=count)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~type+group,nrow = 5) +
      ggplot2::scale_fill_viridis_c(option="magma") +
      ggtitle(names(y_list1)[i]) +
      theme_classic() +
      labs(x = "X (\u03bcm)",y = "Y (\u03bcm)") +
      labs(fill="log-intensity") +
      guides(fill="none")
  })
}

p <- compare_postpred(WW1,WW2,coords)
p
fsave("postpred_example.png")
