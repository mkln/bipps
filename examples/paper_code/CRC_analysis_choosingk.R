devtools::load_all()
library(spatstat)
library(tidyverse)
library(posterior)
library(tidybayes)
library(bayesplot)
library(patchwork)
library(kableExtra)
library(latex2exp)

set.seed(2020)

theme_set(theme_bw(base_size=11, base_family='Times New Roman')+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(angle=45,hjust=1,vjust=1)))

fsave <- \(fname,height=7,width=10) {
  ggsave(paste0(figures_folder,fname),dpi=300, height=height, width=width, units="in")
}
figures_folder <- "examples/data/figures/CRC_choosingk/"

n_thin <- 10
df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots) %>%
  mutate(type = fct_recode(type,"CAFs"="smooth muscle","hybrid E/M"="stroma","TAMs"="CD163+ macros","CTLs"="CD8+ T cells"))
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

# p1 <- plot_y_list(y_list,coords)
# p1
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

        sigma <- unlist(lapply(xl[[k_ix]],\(x) {
          sd(x[ix1,ix2])
        }))

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
ks <- c(2,4,6,8,10)
df1 <- get_df(xl1,out1)
df2 <- get_df(xl2,out2)

lapply(out1,\(o) o$waic)
lapply(out2,\(o) o$waic)

waic_df <- tibble(k=ks,
                  WAIC1=round(unlist(lapply(out1,\(o) o$waic)),0),
                  WAIC2=round(unlist(lapply(out2,\(o) o$waic)),0))

# paper
waic_df %>%
  rename(`WAIC in CLR group`=WAIC1,
         `WAIC in DII group`=WAIC2) %>%
  kable(format = "latex", booktabs = TRUE, digits = 0, align = "c", caption = "The WAIC of the fitted model for various $k$, for both patient groups in the CRC dataset.", label = "waic_crc") %>%
  kable_styling(latex_options = c("striped", "hold_position", "scale_down")) %>%
  add_header_above(c(" " = 1, "WAIC Scores" = 2)) %>%
  column_spec(1, bold = TRUE) %>%
  row_spec(0, bold = TRUE)

# Group 1

# df1 %>%
#   mutate(k=factor(k)) %>%
#   ggplot(aes(hs,ess_bulk,color=k)) +
#   geom_jitter() +
#   facet_wrap(~t1+t2)
# fsave("ess_1.png")

# supplement
df1 %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,rhat,color=k)) +
  geom_jitter() +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  scale_x_continuous(breaks = c(0,1000,2000)) +
  labs(x="Distance (\u03bcm)",y=TeX("$\\hat{R}$"))
fsave("rhat_1.png")

# trace plots
# h_ix <- 1
# trace_df <- lapply(1:length(ks),\(k_ix) {
#   as_draws_df(xl1[[k_ix]][[h_ix]]) %>%
#   pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
#   separate(variable,into = c("type1","type2"),sep=",") %>%
#   mutate(type1 = sub("^x\\[","",type1),
#          type2 = sub("\\]","",type2)) %>%
#   mutate(k=ks[k_ix])
# }) %>%
#   bind_rows()
#
# trace_df %>%
#   mutate(k = factor(k)) %>%
#   mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
#   mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
#   mutate(.chain = factor(.chain)) %>%
#   ggplot(aes(.iteration,value,color=k,group=k)) +
#   geom_line() +
#   geom_hline(yintercept=0,color="black",linetype="dashed",linewidth=0.5) +
#   facet_grid(type1~type2,
#              labeller = label_wrap_gen(width=8)) +
#   # theme_bw() +
#   # theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
#   # theme(axis.title = element_text(size=20),
#   #       axis.text = element_text(size=12),
#   #       strip.text = element_text(size=10)) +
#   labs(x="Draw",y=paste0("Cross-correlation at h=",hs[h_ix]))
# fsave("trace_h0_1.png")


hs <- seq(0,1,0.1)

# lower ones look better, from a convergence standpoint

# supplement
df1 %>%
  # filter(t1 %in% types_intersect[1:3]) %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,E(mu),color=k,fill=k)) +
  geom_line(linewidth=0.5) +
  # geom_ribbon(aes(ymin = lb,ymax=ub)) +
  geom_hline(yintercept = 0,linetype="dotted") +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  scale_x_continuous(breaks = c(0,1000,2000)) +
  scale_y_continuous(breaks = c(-1,0,1)) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")
fsave("xcor_varyingk_group1.png")

# df1 %>%
#   filter(t1 %in% types_intersect[4:10]) %>%
#   mutate(k=factor(k)) %>%
#   ggplot(aes(hs,E(mu),color=k,fill=k)) +
#   geom_line(linewidth=1) +
#   # geom_point() +
#   geom_hline(yintercept = 0,linetype="dotted") +
#   # geom_ribbon(aes(ymin = lb,ymax=ub),alpha=0.5) +
#   facet_wrap(~t1+t2) +
#   theme_minimal() +
#   scale_color_manual(values=as.vector(pals::glasbey()))

# Group 2

# df2 %>%
#   mutate(k=factor(k)) %>%
#   ggplot(aes(hs,ess_bulk,color=k)) +
#   geom_jitter() +
#   facet_wrap(~t1+t2) +
#   theme_minimal()

# supplement
df2 %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,rhat,color=k)) +
  geom_jitter() +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  scale_x_continuous(breaks = c(0,1000,2000)) +
  labs(x="Distance (\u03bcm)",y=TeX("$\\hat{R}$"))
fsave("rhat_2.png")

# trace plots
# h_ix <- 1
# trace_df <- lapply(1:length(ks),\(k_ix) {
#   as_draws_df(xl2[[k_ix]][[h_ix]]) %>%
#     pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
#     separate(variable,into = c("type1","type2"),sep=",") %>%
#     mutate(type1 = sub("^x\\[","",type1),
#            type2 = sub("\\]","",type2)) %>%
#     mutate(k=ks[k_ix])
# }) %>%
#   bind_rows()
#
# trace_df %>%
#   mutate(k = factor(k)) %>%
#   mutate(across(c(type1,type2),~ifelse(.x == "granulocytes","gran.",.x))) %>%
#   mutate(across(c(type1,type2),~ifelse(.x == "vasculature","vasc.",.x))) %>%
#   mutate(.chain = factor(.chain)) %>%
#   ggplot(aes(.iteration,value,color=k,group=k)) +
#   geom_line() +
#   geom_hline(yintercept=0,color="black",linetype="dashed",linewidth=0.5) +
#   facet_grid(type1~type2,
#              labeller = label_wrap_gen(width=8)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
#   theme(axis.title = element_text(size=20),
#         axis.text = element_text(size=12),
#         strip.text = element_text(size=10))

hs <- seq(0,1,0.1)

# lower ones look better, from a convergence standpoint

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
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  scale_x_continuous(breaks = c(0,1000,2000)) +
  scale_y_continuous(breaks = c(-1,0,1)) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  labs(x="Distance (\u03bcm)",y="Cross-correlation")
fsave("xcor_varyingk_group2.png")

# df2 %>%
#   filter(t1 %in% types_intersect[4:10]) %>%
#   mutate(k=factor(k)) %>%
#   ggplot(aes(hs,E(mu),color=k,fill=k)) +
#   geom_line(linewidth=1) +
#   # geom_point() +
#   geom_hline(yintercept = 0,linetype="dotted") +
#   # geom_ribbon(aes(ymin = lb,ymax=ub),alpha=0.5) +
#   facet_wrap(~t1+t2) +
#   theme_minimal() +
#   scale_color_manual(values=as.vector(pals::glasbey()))

# Difference between both groups
group_diff_est <- lapply(1:length(xl1),\(k_ix) {
  lapply(1:length(xl1[[k_ix]]) ,\(i){
    xl1[[k_ix]][[i]] - xl2[[k_ix]][[i]]
  })
})

df_diff <- get_df(group_diff_est,NULL)

# supplement
df_diff %>%
  # filter(t1 %in% types_intersect[1]) %>%
  mutate(k=factor(k)) %>%
  mutate(hs = hs * max_range) %>%
  ggplot(aes(hs,E(mu))) +
  geom_line(aes(color=k)) +
  # geom_ribbon(aes(ymin = lb,ymax=ub,fill=k),alpha=0.5) +
  geom_hline(yintercept = 0,linetype="dotted") +
  # facet_wrap(~t1+t2) +
  facet_grid(t1~t2,
             labeller = labeller(
               t1 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8)),
               t2 = as_labeller(~ recode(., "granulocytes" = "gran.", "vasculature" = "vasc."),
                                label_wrap_gen(width = 8))
             )) +
  scale_x_continuous(breaks = c(0,1000,2000)) +
  # scale_y_continuous(breaks = c(-1,0,1)) +
  scale_color_manual(values=as.vector(pals::glasbey())) +
  # scale_fill_manual(values=as.vector(pals::glasbey())) +
  labs(x="Distance (\u03bcm)",y="Difference in cross-correlation")
fsave("diff_xcor_varyingk.png")
