devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(tidybayes)
library(bayesplot)
library(patchwork)
# library(extrafont)
# font_import()
# loadfonts(device = "win")

set.seed(2020)

theme_bipps <- \() {
  list(
    theme_minimal(),
    theme(
      # text = element_text(family = "Calibri", size = 32, face = "bold"),
      plot.margin = unit(0.5*c(1,1,1,1), "cm")
    )
  )
}

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

df %>%
  ggplot(aes(X,Y,color=Type,fill=Type,shape=Type)) +
  geom_point(size=2.5) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  scale_fill_manual(values=as.vector(pals::glasbey())) +
  scale_shape_manual(values=rep(21:25,times=2,length.out=10)) +
  theme_bipps() +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)")

# binned counts plot
df_raw %>%
  dplyr::group_by(Spot) %>%
  dplyr::mutate(X = X - min(X),
                Y = Y - min(Y)) %>%
  dplyr::ungroup() %>%
  select(X,Y) %>%
  apply(.,2,max) -> max_dim

# dat1 <- df_raw %>%
#   filter(groups == 2) %>%
#   # filter(Spot == "46_B") %>%
#   filter(type %in% types_intersect)

# max_range <- 1919
# x <- dat1$X
# y <- dat1$Y
# types <- dat1$type
# image_ids <- dat1$Spot

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
dplyr::bind_cols(y_list[[1]],coords) %>%
  tidyr::pivot_longer(-c(x,y),names_to = "type",values_to = "Count") %>%
  mutate(x=x*max_range,
         y=y*max_range) %>%
  ggplot2::ggplot(ggplot2::aes(x,y,fill=Count)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~type) +
  scale_fill_viridis_c(option = "inferno",direction=-1,na.value = "transparent") +
  theme_bipps() +
  theme(text = element_text(size=22),
        axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  labs(x="X (\u03BCm)",y="Y (\u03BCm)")

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

out1 <- readRDS("group1_CRC_intercept_lt.rds")
out2 <- readRDS("out2_chains4_CRC_analysis_40k_k2_2e4burn_lt.rds")

lambda <- get_rvars(out1,"lambda")
theta <- get_rvars(out1,"theta")

lambda
hs <- seq(0,1,0.1)

xl1 <- cross_list(out1,hs,thin=40)
xl2 <- cross_list(out2,hs,thin=40)

xl1 <- lapply(xl1,\(xl) {
  rownames(xl) <- colnames(xl) <- types_intersect
  xl
})

xl2 <- lapply(xl2,\(xl) {
  rownames(xl) <- colnames(xl) <- types_intersect
  xl
})

n_types <- length(types_intersect)
n_samples <- ndraws(xl2[[1]])

mcmc_trace(as_draws_df(theta[1,]))

mcmc_trace(as_draws_df(lambda[7:10,]))


# rhat
xl1[[1]] %>%
  summarise_draws() %>%
  separate(variable,into = c("type1","type2"),sep=",") %>%
  mutate(type1 = sub("^\\.\\[","",type1),
         type2 = sub("\\]","",type2)) %>%
  ggplot(aes(type1,type2,fill=rhat)) +
  geom_tile() +
  scico::scale_fill_scico(palette="bam",midpoint=1.1)

# ess_ratio
xl1[[1]] %>%
  summarise_draws(ess=ess_basic) %>%
  separate(variable,into = c("type1","type2"),sep=",") %>%
  mutate(type1 = sub("^\\.\\[","",type1),
         type2 = sub("\\]","",type2)) %>%
  mutate(neff_ratio = ess/n_samples) %>%
  ggplot(aes(type1,type2,fill=neff_ratio)) +
  geom_tile() +
  scico::scale_fill_scico(palette="bam",midpoint=0.1)
# set.seed(2020)
# mask <- matrix(sample(c(TRUE,FALSE),n_types^2,replace=T),nrow=n_types)
# p <- mcmc_trace(as_draws_rvars(xl1[[1]]))
# p

# trace plots
h_ix <- 5
trace_df <- as_draws_df(xl1[[h_ix]]) %>%
  pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
  separate(variable,into = c("type1","type2"),sep=",") %>%
  mutate(type1 = sub("^x\\[","",type1),
         type2 = sub("\\]","",type2))

trace_df %>%
  mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  mutate(.chain = factor(.chain)) %>%
  ggplot(aes(.iteration,value,color=.chain,group=.chain)) +
  geom_line() +
  geom_hline(yintercept=0,color="black",linetype="dashed",linewidth=0.5) +
  facet_grid(type1~type2,
             labeller = label_wrap_gen(width=8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=12),
        strip.text = element_text(size=10)) +
  labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]))

# histogram of draws
trace_df %>%
  mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  mutate(value = ifelse(type1 == type2,NA,value)) %>%
  ggplot(aes(value)) +
  geom_histogram() +
  facet_grid(type1~type2,scales = "free_x")

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
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  labs(x="Type 1",y="Type 2",fill="Correlation") +
  scico::scale_fill_scico(palette = "bam",midpoint = 0) +
  theme_bipps() +
  theme(text=element_text(size=22),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1))

limits <- 2*c(-1,1)
bind_rows(ex1,ex2) %>%
  group_by(type1,type2) %>%
  summarise(diff = value[group == "CLR"] - value[group == "DII"]) %>%
  group_by(type1,type2) %>%
  mutate(lo = quantile(diff,probs = c(0.025)),
         hi = quantile(diff,probs = c(0.975)),
         sig = ifelse(lo < 0 & 0 < hi,FALSE,TRUE)) %>%
  ungroup() %>%
  ggplot(aes(type1,type2,fill=E(diff))) +
  geom_tile() +
  theme_bipps() +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank()) +
  labs(x="Type 1",y="Type 2",fill="Correlation\ndifference") +
  scico::scale_fill_scico(palette = "bam",midpoint = 0,
                          label = function(z) replace(z, c(1, length(z)),
                                                      c(paste0("\u2193 Greater in DII"),
                                                        paste0("\u2191 Greater in CLR"))),
                          breaks = (seq(from=limits[1],
                                        to=limits[2],
                                        length.out=5)),
                          limits = c(limits[1],limits[2]) ) +
  geom_point(aes(shape=ifelse(sig, "dot", "no_dot"),
                                   alpha=ifelse(sig, "dot", "no_dot")),
                      color='red',size=2) +
  scale_shape_manual(values=c(dot=8, no_dot=0), guide="none") +
  scale_alpha_manual(values=c(dot=1, no_dot=0), guide="none") +
  theme(text=element_text(size=22),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1)) +
  theme(plot.margin = unit(0.5*c(1,1,1,1), "cm"))

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

filter_out <- tibble(combo=c("CD163+ macros <--> granulocytes",
                             "CD8+ T cells <--> tumor cells",
                             "B cells <--> granulocytes",
                             "B cells <--> tumor cells",
                             "memory CD4+ T <--> plasma cells",
                             "CD163+ macros <--> tumor cells"))

p1 <- xl1_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  theme_bw() +
  facet_wrap(~combo) +
  # facet_grid(t1~t2,
  #            labeller = label_wrap_gen(width=8)) +
  theme_bipps() +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        title = element_text(size=20),
        plot.margin = unit(0.1*c(1,1,1,1), "cm")) +
  labs(x="",y="Cross-correlation")
  # theme(strip.text = element_text(size=6))

p1

p2 <- xl2_e %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range2) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  theme_bw() +
  facet_wrap(~combo) +
  # facet_grid(t1~t2,
  #            labeller = label_wrap_gen(width=8)) +
  theme_bipps() +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        title = element_text(size=20),
        plot.margin = unit(0.5*c(1,1,1,1), "cm")) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")

p1 / p2 + plot_annotation(tag_levels = 'a')

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

xldiff_e %>%
  # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  mutate(combo = paste0(t1," <--> ",t2)) %>%
  right_join(filter_out) %>%
  mutate(hs = hs * max_range1) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red",linetype="dotted") +
  theme_bipps() +
  # facet_grid(t1~t2,
  #            labeller = label_wrap_gen(width=8)) +
  facet_wrap(~combo) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        title = element_text(size=20)) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation")
