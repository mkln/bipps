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
k_idx <- c(2,3) # corresponds to k=3,4 for PDAC,IPMN
fsave <- \(fname,width=8,height=5) {
  ggsave(paste0(figures_folder,fname),dpi=300, height=height, width=width, units="in")
}
figures_folder <- "examples/data/figures/panc_analysis_paper/"

chains <- 1

dat1 <- readRDS("examples/data/PDAC_POS_DATA.rds") %>%
  rename(X=Cell.X.Position,
         Y=Cell.Y.Position,
         type=Cellname,
         Spot=SlideID) %>%
  filter(type != "CD4")


dat2 <- readRDS("examples/data/IPMN_POS_DATA.rds") %>%
  rename(X=Cell.X.Position,
         Y=Cell.Y.Position,
         type=Cellname,
         Spot=SlideID) %>%
  filter(type != "CD4")

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

dat1 %>%
  count(Spot)

df <- dat1 %>%
  filter(Spot == "PATID_64_SLIDE_1")
  # filter(Spot == "PATID_12_SLIDE_1")
  # filter(Spot == "PATID_538_SLIDE_1")

# pp plot
df %>%
  ggplot(aes(X,Y,color=type,fill=type,shape=type)) +
  geom_point(size=1) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  scale_fill_manual(values=as.vector(pals::glasbey())) +
  scale_shape_manual(values=rep(21:25,times=2,length.out=10)) +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)",color="Type",shape="Type",fill="Type")
fsave("pdac_pat.png")

# binned counts plot
bind_rows(dat1,dat2) %>%
  dplyr::group_by(Spot) %>%
  dplyr::mutate(X = X - min(X),
                Y = Y - min(Y)) %>%
  dplyr::ungroup() %>%
  select(X,Y) %>%
  apply(.,2,max) -> max_dim

pix_dim <- 70
nx <- ceiling(max_dim[1]/pix_dim)
ny <- ceiling(max_dim[2]/pix_dim)

out <- pixellate_grid(df$X,df$Y,df$type,df$Spot,nx,ny)

y_list <- out$y_list
coords <- out$coords

# sum(is.na(y_list[[1]][,1])) / nrow(coords)

max_range <- max(max_dim)
dplyr::bind_cols(y_list[[1]],coords) %>%
  tidyr::pivot_longer(-c(x,y),names_to = "type",values_to = "Count") %>%
  mutate(x=x*max_range,
         y=y*max_range) %>%
  ggplot2::ggplot(ggplot2::aes(x,y,fill=Count)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~type) +
  scale_fill_viridis_c(option = "inferno",direction=-1,na.value = "transparent",breaks = c(0,5,10)) +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)")
fsave("binned_count_PDAC.png")

out2 <- readRDS("examples/data/IPMN_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")
out1 <- readRDS("examples/data/PDAC_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")

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

xl1 <- xl1[[k_idx[1]]]
xl2 <- xl2[[k_idx[2]]]


# zero distance expected spatial cross-cor plots
ex1 <- xl1[[1]] %>%
  as_tibble(rownames="type1") %>%
  pivot_longer(-type1,names_to = "type2") %>%
  mutate(group = "PDAC")

ex2 <- xl2[[1]] %>%
  as_tibble(rownames="type1") %>%
  pivot_longer(-type1,names_to = "type2") %>%
  mutate(group = "IPMN")

bind_rows(ex1,ex2) %>%
  ggplot(aes(type1,type2,fill=E(value))) +
  geom_tile() +
  facet_wrap(~group) +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  labs(x="Type 1",y="Type 2",fill="Correlation") +
  scico::scale_fill_scico(palette = "bam",midpoint = 0) +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))
fsave("exp_zero_cor.png")

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
  ungroup() %>%
  mutate(group = "PDAC") -> xl1_e

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
  ungroup() %>%
  mutate(group = "IPMN") -> xl2_e

p1 <- bind_rows(xl1_e,xl2_e) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  mutate(hs = hs * max_range) %>%
  filter(t1 == t2) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub,group=group),fill = "grey70") +
  geom_line(aes(color=group)) +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  facet_wrap(~t2,nrow=1) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  labs(x="Distance (\u03bcm)",y="Marginal correlation") +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())
p1

p2 <- bind_rows(xl1_e,xl2_e) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  mutate(hs = hs * max_range) %>%
  filter(t1 != t2) %>%
  filter(t1 != "PDL1_CD8",
         t2 != "PDL1_CD8") %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub,group=group),fill = "grey70") +
  geom_line(aes(color=group)) +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  facet_wrap(~combo,nrow=1) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  guides(color="none") +
  theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

p3 <- bind_rows(xl1_e,xl2_e) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  mutate(hs = hs * max_range) %>%
  filter(t1 == "PDL1_CD8" |
         t2 == "PDL1_CD8") %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub,group=group),fill = "grey70") +
  geom_line(aes(color=group)) +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  facet_wrap(~combo,nrow=1) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation") +
  guides(color="none")

p1/p2/p3 + plot_annotation(tag_levels = "a")


fsave("xcor_comparison.png",width=9,height=6)

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

# paper
xldiff_e %>%
  # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  # right_join(filter_out) %>%
  mutate(hs = hs * max_range) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  facet_grid(t1~t2,
             labeller = label_wrap_gen(width=8)) +
  # facet_wrap(~combo) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=14),
  #       strip.text = element_text(size=16),
  #       title = element_text(size=20)) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation")
fsave("diff_cor_panc.png")
