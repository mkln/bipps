rm(list=ls())
devtools::load_all()
library(magrittr)
library(spatstat)
library(pals)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(2020)

df_raw <- readr::read_csv("examples/data/CRC_cleaned.csv") %>%
  dplyr::mutate(type = as.factor(type)) %>%
  dplyr::rename(Spot = spots)

spot <- "59_A"

df = df_raw %>%
  dplyr::filter(Spot == spot) %>%
  dplyr::mutate(type = forcats::fct_lump_n(type,3)) %>%
  dplyr::filter(type != "Other") %>%
  droplevels()

window <- owin(xrange = c(min(df$X),max(df$X)), yrange = c(min(df$Y),max(df$Y)))
pp <- ppp(df$X,df$Y,marks=df$type,window = window)

pp %>%
  as.data.frame() %>%
  ggplot(aes(x,y,color=marks,shape=marks)) +
  geom_point() +
  scale_color_manual(values = as.vector(glasbey(nlevels(df$type)))) +
  scale_shape_manual(values = 1:nlevels(df$type))


midpoint <- function(x) {
  x[1] + (x[2]-x[1])/2
}

extract_midpoints <- Vectorize(function(x) {
  midpoint(as.numeric(stringr::str_split(gsub("\\[|)|\\]","",x),",")[[1]]))
})

nx <- 20

qd <- quadratcount(split(pp),nx=nx)

plot(qd)

qd <- qd %>%
  as.data.frame() %>%
  rename(x=2,y=1) %>%
  select(-contains(".x"),-contains(".y")) %>%
  rename_with(~gsub(".Freq","",.x),-c(x,y)) %>%
  mutate(across(c(x,y),extract_midpoints))

coords <- qd %>%
  select(x,y) %>%
  as.matrix()

YY <- qd %>%
  select(-c(x,y)) %>%
  as.matrix()

XX <- matrix(0,nrow=nrow(YY),ncol=ncol(YY))
# XX <- NULL


mcmc_keep <- 100
mcmc_burn <- 500
mcmc_thin <- 1

q <- ncol(YY)
k <- 3
Lambda <- matrix(0, q, k)

set.seed(1)
mesh_total_time <- system.time({
  meshout <- bipps(YY,XX, coords, k = k,
                   family="poisson",
                   block_size=5,
                   n_samples = mcmc_keep, n_burn = mcmc_burn, n_thin = mcmc_thin,
                   n_threads = 2,
                   starting=list(lambda = Lambda, phi=1),
                   prior = list(phi=c(1, 200), nu=c(.5, .5)),
                   settings = list(adapting=T, saving=F, ps=T, hmc=0),
                   verbose=10,
                   debug=list(sample_beta=F, sample_tausq=F,
                              sample_theta=F, sample_w=T, sample_lambda=T,
                              verbose=F, debug=T)
  )})


plot_cube <- function(cube_mcmc, q, k, name="Parameter"){
  par(mar=c(2.5,2,1,1), mfrow=c(q,k))
  for(i in 1:q){
    for(j in 1:k){
      cube_mcmc[i, j,] %>% plot(type='l', main="{name} {i}, {j}" %>% glue::glue())
    }
  }
}

# chain plots
plot_cube(meshout$theta_mcmc, 1, k, "theta")
plot_cube(meshout$lambda_mcmc, q, k, "Lambda")
# plot_cube(meshout$beta_mcmc, p, q, "Beta")

# posterior means
# meshout$tausq_mcmc %>% apply(1, mean)
(lambda <- meshout$lambda_mcmc %>% apply(1:2, mean))
# meshout$beta_mcmc %>% apply(1:2, mean)
meshout$theta_mcmc %>% apply(1:2, mean)

# what are theta, tausq in this example? Can we just remove covariates if we don't have them?
# why is theta end up fixed to 1? and have a zero% acceptance rate? why is tausq fixed to 1?
# where does phi show up?
# why do covariates have to be column size equal to number of outcomes?

# process means
wmesh <- data.frame(meshout$w_mcmc %>% summary_list_mean())
colnames(wmesh) <- paste0("wmesh_", 1:q)
# predictions
ymesh <- data.frame(meshout$yhat_mcmc %>% summary_list_mean())
colnames(ymesh) <- paste0("pred_", colnames(YY))

# observed
YY_obs <- YY
colnames(YY_obs) <- paste0("obs_",colnames(YY))
mesh_df <-
  meshout$coordsdata %>%
  cbind(ymesh) %>%
  # cbind(wmesh) %>%
  cbind(YY_obs)

res <- mesh_df %>%
  pivot_longer(-c(x,y),names_pattern = "(pred|obs)_(.*)$",names_to=c(".value","type"))


res %>%
    pivot_longer(c(pred,obs)) %>%
    ggplot(aes(x,y, fill=value)) +
    geom_tile() +
    facet_wrap(name ~ type, ncol= q) +
    scico::scale_fill_scico(palette = "batlow")

res %>%
  group_by(type) %>%
  summarise(rmse=sqrt(mean((pred-obs)^2)))

(lower <- meshout$lambda_mcmc %>%
  apply(.,1:2,\(x) quantile(x,probs=c(0.025))))

(upper <- meshout$lambda_mcmc %>%
  apply(.,1:2,\(x) quantile(x,probs=c(0.975))))
