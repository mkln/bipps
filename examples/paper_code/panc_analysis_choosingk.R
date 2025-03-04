devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(patchwork)
library(kableExtra)
library(ggdist)
library(latex2exp)

set.seed(2020)

theme_set(theme_bw(base_size=11, base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle=45,hjust=1,vjust=1)))

fsave <- \(fname) {
  ggsave(paste0(figures_folder,fname),dpi=300, height=5, width=8, units="in")
}
figures_folder <- "examples/data/figures/panc_choosingk/"

n_samples <- 1000
n_burnin <- 20000
n_thin <- 10
n_threads <- 8
block_size <- 50
ks <- c(2,3,4,5)
starting <- list(phi = 5)
prior <- list(phi = c(0.1,10))

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

max_range <- max(max_dim)


out1 <- readRDS("examples/data/PDAC_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")
out2 <- readRDS("examples/data/IPMN_varyingk_nburn20k_nthin10_nsamp1e3_chain1_lt.rds")

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

unique_combinations_with_self <- function(data) {
  # Generate all pairs including self-pairings
  all_pairs <- expand.grid(data, data)

  # Sort the pairs and remove duplicates
  all_pairs <- t(apply(all_pairs, 1, sort))
  unique_pairs <- unique(all_pairs)

  return(unique_pairs)
}

get_df <- \(xl,out) {

  lapply(1:length(xl),\(k_ix) {
    unique_combinations_with_self(types_intersect) %>%
      as.data.frame() %>%
      as_tibble() %>%
      magrittr::set_colnames(c("t1","t2")) %>%
      group_by(t1,t2) %>%
      group_modify(~{
        ix1 <- which(types_intersect == .y$t1)
        ix2 <- which(types_intersect == .y$t2)

        mu <- as.vector(do.call(rbind,lapply(xl[[k_ix]],\(x) {
          x[ix1,ix2]
        })))

        ess_bulk <- as.vector(do.call(rbind,lapply(xl[[k_ix]],\(x) {
          posterior::ess_bulk(x[ix1,ix2])
        })))

        ess_tail <- as.vector(do.call(rbind,lapply(xl[[k_ix]],\(x) {
          posterior::ess_tail(x[ix1,ix2])
        })))

        rhats <- as.vector(do.call(rbind,lapply(xl[[k_ix]],\(x) {
          posterior::rhat(x[ix1,ix2])
        })))

        lb <- unlist(lapply(xl[[k_ix]],\(x) {
          quantile(x[ix1,ix2],probs = 0.025)
        }))

        ub <- unlist(lapply(xl[[k_ix]],\(x) {
          quantile(x[ix1,ix2],probs = 0.975)
        }))
        tibble(mu=mu,lb=lb,ub=ub,hs=hs,
               ess_bulk=ess_bulk,ess_tail=ess_tail,rhat=rhats)
      }) %>%
      ungroup() %>%
      mutate(k=ks[k_ix]) -> df

    if(!is.null(out)) {
      df$waic <- out[[k_ix]]$waic
    }

    return(df)
  }) %>%
    bind_rows() -> df

}

df1 <- get_df(xl1,out1)
df2 <- get_df(xl2,out2)

lapply(out1,\(o) o$waic)
lapply(out2,\(o) o$waic)

waic_df <- tibble(k=ks,
                  WAIC1=round(unlist(lapply(out1,\(o) o$waic)),0),
                  WAIC2=round(unlist(lapply(out2,\(o) o$waic)),0))

# paper
waic_df %>%
  rename(`WAIC in PDAC group`=WAIC1,
         `WAIC in IPMN group`=WAIC2) %>%
  kable(format = "latex", booktabs = TRUE, digits = 0, align = "c", caption = "", label = "waic_panc") %>%
  kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
  add_header_above(c(" " = 1, "WAIC Scores" = 2)) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE)

df1 %>%
  mutate(hs = hs * max_range) %>%
  mutate(k=factor(k)) %>%
  ggplot(aes(hs,rhat,color=k)) +
  geom_jitter() +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2) +
  labs(x="Distance (\u03bcm)",y=TeX("$\\hat{R}$"))
fsave("rhat_1.png")

# trace plots
h_ix <- 1
trace_df <- lapply(1:length(ks),\(k_ix) {
  as_draws_df(xl1[[k_ix]][[h_ix]]) %>%
    pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
    separate(variable,into = c("type1","type2"),sep=",") %>%
    mutate(type1 = sub("^x\\[","",type1),
           type2 = sub("\\]","",type2)) %>%
    mutate(k=ks[k_ix])
}) %>%
  bind_rows()

trace_df %>%
  # filter(k == 6) %>%
  mutate(k = factor(k)) %>%
  # mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  mutate(.chain = factor(.chain)) %>%
  ggplot(aes(.iteration,value,color=k,group=k)) +
  geom_line() +
  geom_hline(yintercept=0,color="black",linetype="dashed",linewidth=0.5) +
  facet_grid(type1~type2,
             labeller = label_wrap_gen(width=8)) +
  # theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))
  # theme(axis.title = element_text(size=20),
  #       axis.text = element_text(size=12),
  #       strip.text = element_text(size=10))

df1 %>%
  # filter(t1 %in% types_intersect[1:3]) %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,E(mu),color=k,fill=k)) +
  geom_line(linewidth=0.5) +
  # geom_ribbon(aes(ymin = lb,ymax=ub)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2) +
  # theme_minimal() +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")
fsave("xcor_varyingk_group1.png")

df2 %>%
  mutate(hs = hs * max_range) %>%
  mutate(k=factor(k)) %>%
  ggplot(aes(hs,rhat,color=k)) +
  geom_jitter() +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2) +
  labs(x="Distance (\u03bcm)",y=TeX("$\\hat{R}$"))
fsave("rhat_2.png")

df2 %>%
  group_by(t1,t2,k) %>%
  summarise(mn_rhat = mean(rhat)) %>%
  ggplot(aes(k,y=mn_rhat)) +
  geom_point() +
  facet_grid(t1~t2)

df2 %>%
  group_by(k) %>%
  summarise(mn_rhat = mean(rhat)) %>%
  ggplot(aes(k,y=mn_rhat)) +
  geom_point()


h_ix <- 5
trace_df <- lapply(1:length(ks),\(k_ix) {
  as_draws_df(xl2[[k_ix]][[h_ix]]) %>%
    pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
    separate(variable,into = c("type1","type2"),sep=",") %>%
    mutate(type1 = sub("^x\\[","",type1),
           type2 = sub("\\]","",type2)) %>%
    mutate(k=ks[k_ix])
}) %>%
  bind_rows()

trace_df %>%
  filter(k == 5) %>%
  mutate(k = factor(k)) %>%
  # mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
  # mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
  mutate(.chain = factor(.chain)) %>%
  ggplot(aes(.iteration,value,color=k,group=k)) +
  geom_line() +
  geom_hline(yintercept=0,color="black",linetype="dashed",linewidth=0.5) +
  facet_grid(type1~type2,
             labeller = label_wrap_gen(width=8)) +
  # theme_bw() +
  theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1))

# supplement
df2 %>%
  # filter(t1 %in% types_intersect[1:3]) %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,E(mu),color=k,fill=k)) +
  geom_line(linewidth=0.5) +
  # geom_ribbon(aes(ymin = lb,ymax=ub)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2) +
  # theme_minimal() +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")
fsave("xcor_varyingk_group2.png")
