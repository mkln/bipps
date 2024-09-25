rm(list = ls())
devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(tidybayes)
library(bayesplot)

set.seed(2020)

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
  theme_bw() +
  theme(text=element_text(size=32)) +
  theme(plot.margin = unit(0.5*c(1,1,1,1), "cm"))

# binned counts plot
df_raw %>%
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
out <- create_y_list(df$X,df$Y,df$Type,df$Spot,nx,ny)
y_list <- out$y_list
coords <- out$coords

dplyr::bind_cols(y_list[[1]],coords) %>%
  tidyr::pivot_longer(-c(x,y),names_to = "type",values_to = "Count") %>%
  ggplot2::ggplot(ggplot2::aes(x,y,fill=Count)) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~type) +
  scico::scale_fill_scico() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),panel.grid.major=element_blank(),
        axis.text = element_blank()) +
  labs(x="",y="") +
  theme(text=element_text(size=22)) +
  theme(plot.margin = unit(0.5*c(1,1,1,1), "cm"))

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

out1 <- readRDS("out1_chains_CRC_analysis_v2.rds")
out2 <- readRDS("out2_chains_CRC_analysis.rds")

lambda <- get_rvars(out1,"lambda",thin=40)
theta <- get_rvars(out1,"theta",thin=40)

lambda
hs <- seq(0,100,10)

xl1 <- cross_list(out1,hs,thin=40)
xl2 <- cross_list(out2,hs)

xl1 <- lapply(xl1,\(xl) {
  rownames(xl) <- colnames(xl) <- types_intersect
  xl
})

xl2 <- lapply(xl2,\(xl) {
  rownames(xl) <- colnames(xl) <- types_intersect
  xl
})

n_types <- length(types_intersect)
n_samples <- ndraws(xl1[[1]])

summarise_draws(lambda)

summarise_draws(xl1[[2]])

mcmc_trace(as_draws_df(theta[1,]))

mcmc_trace(as_draws_df(lambda[7:8,]))


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
h_ix <- 3
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
  theme(text=element_text(size=22),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1)) +
  theme(plot.margin = unit(0.5*c(1,1,1,1), "cm"))

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
  theme_minimal() +
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


xl1_e %>%
  # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
  # mutate(ub = mu+sigma,
  #        lb = mu-sigma) %>%
  mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red") +
  theme_bw() +
  facet_grid(t1~t2,
             labeller = label_wrap_gen(width=8)) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=12),
        strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")
  # theme(strip.text = element_text(size=6))

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
  mutate(across(c(t1,t2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  mutate(across(c(t1,t2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  ggplot(aes(hs,mu)) +
  geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
  geom_line() +
  geom_hline(yintercept = 0,color="red") +
  theme_bw() +
  facet_grid(t1~t2,
             labeller = label_wrap_gen(width=8)) +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=12),
        strip.text = element_text(size=10)) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation")
