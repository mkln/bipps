devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(kableExtra)

set.seed(2020)

theme_set(theme_bw(base_size=11, base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank()))
n_thin <- 10
k_idx <- 3 # corresponds to k=6
fsave <- \(fname,height=5,width=8) {
  ggsave(paste0(figures_folder,fname),dpi=300, height=height, width=width, units="in")
}
figures_folder <- "examples/data/figures/CRC_analysis_paper/"

df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))
# mutate(Spot = factor(Spot))

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

# pp plot
df <- df_raw %>%
  filter(type %in% types_intersect) %>%
  filter(Spot == "59_A") %>%
  # filter(Spot == "67_B") %>%
  # filter(Spot == "57_A") %>%
  # filter(Spot == "53_B") %>%
  # filter(Spot == "54_B") %>%
  rename(Type=type)

df %>%
  ggplot(aes(X,Y,color=Type,fill=Type,shape=Type)) +
  geom_point(size=2.5) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  scale_fill_manual(values=as.vector(pals::glasbey())) +
  scale_shape_manual(values=rep(21:25,times=2,length.out=10)) +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)")
fsave("crc_pat.png")

# binned counts plot
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

# sum(!is.na(y_list[[1]][,1])) / nrow(coords)
# 53%, one cell type
# 47%, one cell type
# 51%, two cell types
# 63%, 4 cell types

max_range <- max(max_dim)
dplyr::bind_cols(y_list[[1]],coords) %>%
  tidyr::pivot_longer(-c(x,y),names_to = "type",values_to = "Count") %>%
  mutate(x=x*max_range,
         y=y*max_range) %>%
  ggplot2::ggplot(ggplot2::aes(x,y,fill=Count)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~type,nrow=2) +
  scale_fill_viridis_c(option = "inferno",direction=-1,na.value = "transparent") +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)")
fsave("binned_count_CRC.png")

out1 <- readRDS("examples/data/group1_CRC_intercept_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")
out2 <- readRDS("examples/data/group2_CRC_intercept_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")

# so it is between 8 and 10 factors, according to WAIC.

hs <- seq(0,1,0.1)

xl1 <- lapply(out1,\(o) cross_list(list(o),hs,thin=n_thin))
xl2 <- lapply(out2,\(o) cross_list(list(o),hs,thin=n_thin))

xl1 <- lapply(xl1,\(xl) {
  lapply(xl,\(x) {
    rownames(x) <- colnames(x) <- types_intersect
    x
  })
})

xl2 <- lapply(xl2,\(xl) {
  lapply(xl,\(x) {
    rownames(x) <- colnames(x) <- types_intersect
    x
  })
})

xl1 <- xl1[[k_idx]]
xl2 <- xl2[[k_idx]]


# zero distance expected spatial cross-cor plots
ex1 <- xl1[[1]] %>%
  as_tibble(rownames="type1") %>%
  pivot_longer(-type1,names_to = "type2") %>%
  mutate(group = "CLR")

ex2 <- xl2[[1]] %>%
  as_tibble(rownames="type1") %>%
  pivot_longer(-type1,names_to = "type2") %>%
  mutate(group = "DII")

bind_rows(ex1,ex2) %>%
  ggplot(aes(type1,type2,fill=E(value))) +
  geom_tile() +
  facet_wrap(~group) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  labs(x="Type 1",y="Type 2",fill="Correlation") +
  scico::scale_fill_scico(palette = "bam",midpoint = 0) +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
fsave("exp_zero_cor_CRC.png")


unique_combinations_with_self <- function(data) {
  # Generate all pairs including self-pairings
  all_pairs <- expand.grid(data, data)

  # Sort the pairs and remove duplicates
  all_pairs <- t(apply(all_pairs, 1, sort))
  unique_pairs <- unique(all_pairs)

  return(unique_pairs)
}

# spatial cross-correlation over distance plots
unique_combinations_with_self(types_intersect) %>%
  # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types_intersect == .y$t1)
    ix2 <- which(types_intersect == .y$t2)

    mu <- unlist(lapply(xl1,\(x) {
      E(x[ix1,ix2])
    }))

    int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

    sigma <- unlist(lapply(xl1,\(x) {
      sd(x[ix1,ix2])
    }))

    lb <- unlist(lapply(xl1,\(x) {
      quantile(x[ix1,ix2],probs = 0.025)
    }))

    ub <- unlist(lapply(xl1,\(x) {
      quantile(x[ix1,ix2],probs = 0.975)
    }))
    tibble(mu=mu,lb=lb,ub=ub,hs=hs,auc=int_val)
  }) %>%
  ungroup() -> xl1_e

# spatial cross-correlation over distance plots
unique_combinations_with_self(types_intersect) %>%
  # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("t1","t2")) %>%
  group_by(t1,t2) %>%
  group_modify(~{
    ix1 <- which(types_intersect == .y$t1)
    ix2 <- which(types_intersect == .y$t2)

    mu <- unlist(lapply(xl2,\(x) {
      E(x[ix1,ix2])
    }))

    int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

    sigma <- unlist(lapply(xl2,\(x) {
      sd(x[ix1,ix2])
    }))

    lb <- unlist(lapply(xl2,\(x) {
      quantile(x[ix1,ix2],probs = 0.025)
    }))

    ub <- unlist(lapply(xl2,\(x) {
      quantile(x[ix1,ix2],probs = 0.975)
    }))
    tibble(mu=mu,lb=lb,ub=ub,hs=hs,auc=int_val)
  }) %>%
  ungroup() -> xl2_e


filter_out <- tibble(combo=c("hybrid E/M <--> TAMs",
                             "CTLs <--> tumor cells",
                             "B cells <--> granulocytes",
                             "B cells <--> tumor cells",
                             "memory CD4+ T <--> plasma cells",
                             "CAFs <--> CTLs"))

dd1 <- xl1_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range) %>%
  mutate(group = "CLR")

dd2 <- xl2_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range) %>%
  mutate(group = "DII")


bind_rows(dd1,dd2) %>%
  rename(Group = group) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub,fill = Group),alpha=0.2) +
  geom_line(aes(color=Group)) +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  facet_wrap(~combo) +
  # facet_grid(t1~t2,
  #            labeller = label_wrap_gen(width=8)) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=14),
  #       strip.text = element_text(size=16),
  #       title = element_text(size=20),
  #       plot.margin = unit(0.5*c(1,1,1,1), "cm")) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")

# paper
fsave("xcor_CRC.png")

# supplement - full version
dd1 <- xl1_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  mutate(hs = hs * max_range) %>%
  mutate(group = "CLR")

dd2 <- xl2_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  mutate(hs = hs * max_range) %>%
  mutate(group = "DII")


bind_rows(dd1,dd2) %>%
  rename(Group = group) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub,fill = Group),alpha=0.2) +
  geom_line(aes(color=Group)) +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  # facet_wrap(~combo) +
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=14),
  #       strip.text = element_text(size=16),
  #       title = element_text(size=20),
  #       plot.margin = unit(0.5*c(1,1,1,1), "cm")) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  scale_x_continuous(breaks = c(0,1000,2000))
fsave("xcor_CRC_full.png",height = 7,width = 10)

# spatial cross-correlation difference over distance plots
xl_diff <- lapply(1:length(hs),\(i) {
  x1 <- xl1[[i]]
  x2 <- xl2[[i]]

  x1 - x2
})

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

# supplement
xldiff_e %>%
  # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  # facet_grid(t1~t2,
  #            labeller = label_wrap_gen(width=8)) +
  facet_wrap(~combo) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=14),
  #       strip.text = element_text(size=16),
  #       title = element_text(size=20)) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation")
fsave("diff_cor_CRC.png")

# supplement - full version
