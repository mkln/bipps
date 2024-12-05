# spatially stationary W
# (also try with nonstationary)
# vary -> number of outcomes, spatial factors and images. lambda and theta.

devtools::load_all()
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(spatstat)
library(fields)

set.seed(2020)

# mcmc settings
n_samples <- 500
n_burnin <- 20000
n_thin <- 20
n_threads <- 4
block_size <- 50
Theta <- 0.7 # scaled version
starting <- list(phi = Theta)
prior <- list(phi = c(0.1,10))
chains <- 1
do_plots <- FALSE
save_file <- "out_sim3.rds"
save_file_lt <- "out_sim3_lt.rds"

# simulation settings
nx <- 30
ny <- 30
num_images <- 20
# theta
x_max <- 1919
y_max <- 1439
inv_theta <- 1 / Theta * max(x_max,y_max)
sigmasq <- 1

scaling <- 20
mu <- 0
k <- 2
q <- 2
p <- 1

grid<- list( x= seq( 0,x_max*scaling,,nx*scaling), y= seq(0,y_max*scaling,,ny*scaling))
obj<- circulantEmbeddingSetup( grid, Covariance="Exponential", theta=inv_theta)

# subset V
idx_x <- which(grid[[1]] < x_max)
idx_y <- which(grid[[2]] < y_max)

gridx <- grid[[1]][idx_x]
gridy <- grid[[2]][idx_y]

set.seed(2025)
V <- lapply(1:num_images,\(i) {
  lapply(1:k,\(i) {
    sim<-  sqrt(sigmasq)*circulantEmbedding( obj)
    subW <- sim[idx_y,idx_x]

    print(range(subW))
    subW
  })
})

image.plot(gridx,gridy,V[[1]][[1]])


VV <- lapply(1:num_images,\(i) {
  lapply(V[[i]],\(w) c(w)) %>%
    do.call(cbind,.)
})


# factor loadings
set.seed(2025)
Lambda <- matrix(0, q, k)
diag(Lambda) <- runif(k, 0.5, 1)
Lambda[lower.tri(Lambda)] <- runif(sum(lower.tri(Lambda)), -0.7, 0.7)

# nuggets
# tau.sq <- rep(.01, q)
# TTsq <- matrix(1, nrow=n) %x% matrix(tau.sq, ncol=length(tau.sq))
# # measurement errors
# # double check these
# EE <- ( rnorm(n*length(tau.sq)) %>% matrix(ncol=length(tau.sq)) ) * TTsq^.5

# XX <- lapply(1:num_images,\(j) {
#   1:p %>% lapply(function(i) rnorm(nrow(VV[[i]]), 0, .1)) %>% do.call(cbind, .)
# })
# Beta <- matrix(rnorm(p*q), ncol=q) * 0

WW <- lapply(1:num_images,\(i) {
  mat <- VV[[i]] %*% t(Lambda) + mu
  print(range(exp(mat)))
  mat
})

# expand.grid(x=gridx,y=gridy) %>%
#   mutate(value=WW[[1]][,1]) %>%
#   ggplot(aes(x,y,fill=value)) +
#   geom_tile()
sz_x <- length(gridx)
sz_y <- length(gridy)
W <- lapply(1:num_images,\(i) {
  lapply(1:q,\(j) {
    mat <- matrix(WW[[i]][,j],nrow=sz_y,ncol=sz_x)
    print(range(exp(mat)))
    mat
  })
})

image.plot(gridx,gridy,W[[1]][[2]])
#
# # make linear predictor
# lin_pred <- lapply(1:num_images,\(i) {
#   lapply(1:q,\(j) {
#     mat <- W[[i]][[j]] + mu #+ matrix(XX[[i]][,j] %*% Beta,nrow=sz_y,ncol=sz_x)
#     print(range(mat))
#     mat
#   })
# })

y_list <- lapply(WW,\(ww) {
  matrix(rpois(nrow(ww)*ncol(ww),exp(ww)),nrow=nrow(ww),ncol=ncol(ww))
})

coords <- expand.grid(x=gridx,y=gridy)

coords <- coords / max(x_max,y_max)


# make point patterns
# window <- owin(c(0,gridx[length(gridx)]),c(0,gridy[length(gridy)]))
#
# imgs <- lapply(1:num_images,\(i) {
#   lapply(1:q,\(j) {
#     img <- as.im(t(exp(lin_pred[[i]][[j]])),W=window)
#     print(mean(img)*x_max*y_max)
#     img
#   })
# })

# img <- imgs[[1]][[1]]
# mean(img)*x_max*y_max
#
# pats <- lapply(1:num_images,\(i) {
#   pp <- rmpoispp(imgs[[i]])
# })
#
# pp_df <- lapply(1:num_images,\(i) {
#   pats[[i]] %>%
#     as.data.frame() %>%
#     rename(type = marks) %>%
#     mutate(spot = paste0("spot_",i))
# }) %>%
#   bind_rows()

# pp_df %>%
#   ggplot(aes(x,y,color=type)) +
#   geom_point() +
#   theme_classic() +
#   facet_wrap(~spot,ncol=1)

# pixellate the point patterns
# out <- pixellate_grid(pp_df$x,pp_df$y,pp_df$type,pp_df$spot,28,21)
#
# y_list <- out$y_list
# coords <- out$coords
x_list <- lapply(y_list,\(yy) {
  matrix(0,nrow = nrow(yy),ncol = 1)
})

if(do_plots) {
  p1 <- plot_y_list(y_list,coords)
  p1
}

# run bipps
out <- lapply(1:chains,\(i) multi_bipps(y_list,
                                         x_list,
                                         coords,
                                         k = k,
                                         family = "poisson",
                                         block_size = block_size,
                                         n_samples = n_samples, n_burn = n_burnin, n_thin = n_thin,
                                         n_threads = n_threads,
                                         starting = starting,
                                         prior = prior,
                                         settings = list(adapting = T, saving = T, ps = T,low_mem=T),
                                         verbose = 10,
                                         debug = list(
                                           sample_beta = T, sample_tausq = F,
                                           sample_theta = T, sample_w = T, sample_lambda = T,
                                           verbose = F, debug = F
                                         ),
                                         just_preprocess = F))
saveRDS(out,save_file)

out_lt <- lapply(out,\(o) list(theta_mcmc=o$theta_mcmc,lambda_mcmc=o$lambda_mcmc))

saveRDS(out_lt,save_file_lt)

if(do_plots) {
  out <- readRDS("out_sim3_lt.rds")

  lambda <- get_rvars(out,"lambda",thin=n_thin)
  lambda
  mcmc_trace(as_draws_df(lambda[,2]))
  theta <- get_rvars(out,"theta",thin=n_thin)
  theta
  mcmc_trace(as_draws_df(theta[1,]))

  hs <- seq(0,1,0.1)
  xl <- cross_list(out,hs,thin=n_thin)
  out_actual <- list(theta_mcmc=matrix(c(rep(Theta,k),rep(0,k)),nrow = 2,ncol=k,byrow=TRUE),
                     lambda_mcmc=Lambda)
  out_actual <- list(lapply(out_actual,\(o) {
    dim(o) <- c(dim(o),1)
    o
  }))

  h_ix <- 1
  trace_df <- as_draws_df(xl[[h_ix]]) %>%
    pivot_longer(-c(".chain",".iteration",".draw"),names_to = "variable") %>%
    separate(variable,into = c("type1","type2"),sep=",") %>%
    mutate(type1 = sub("^x\\[","",type1),
           type2 = sub("\\]","",type2))

  trace_df %>%
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

  xl_actual <- cross_list(out_actual,hs,thin=n_thin)
  xl_actual <- lapply(xl_actual,\(x) E(x))

  unique_combinations_with_self <- function(data) {
    # Generate all pairs including self-pairings
    all_pairs <- expand.grid(data, data)

    # Sort the pairs and remove duplicates
    all_pairs <- t(apply(all_pairs, 1, sort))
    unique_pairs <- unique(all_pairs)

    return(unique_pairs)
  }

  types <- 1:q
  # spatial cross-correlation over distance plots
  unique_combinations_with_self(types) %>%
    # expand_grid(type1 = types_intersect,type2 = types_intersect) %>%
    as.data.frame() %>%
    as_tibble() %>%
    magrittr::set_colnames(c("t1","t2")) %>%
    group_by(t1,t2) %>%
    group_modify(~{
      ix1 <- which(types == .y$t1)
      ix2 <- which(types == .y$t2)

      mu <- unlist(lapply(xl,\(x) {
        E(x[ix1,ix2])
      }))

      mu_actual <- unlist(lapply(xl_actual,\(x) {
        x[ix1,ix2]
      }))

      int_val <- integrate(approxfun(hs,abs(mu)),hs[1],hs[length(hs)])$value

      sigma <- unlist(lapply(xl,\(x) {
        sd(x[ix1,ix2])
      }))

      lb <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.025)
      }))

      ub <- unlist(lapply(xl,\(x) {
        quantile(x[ix1,ix2],probs = 0.975)
      }))
      tibble(mu=mu,mu_actual,lb=lb,ub=ub,hs=hs,auc=int_val)
    }) %>%
    ungroup() -> xl_e


  xl_e %>%
    mutate(hs = hs*x_max) %>%
    pivot_longer(c(mu,mu_actual)) %>%
    # filter(auc %in% sort(unique(auc),decreasing = TRUE)[1:10]) %>%
    # mutate(ub = mu+sigma,
    #        lb = mu-sigma) %>%
    ggplot() +
    geom_ribbon(aes(hs,ymin = lb,ymax=ub),fill = "grey70") +
    geom_line(aes(hs,value,color=name)) +
    geom_hline(yintercept = 0,color="red",alpha=0.5,linetype="dotted") +
    theme_minimal() +
    facet_grid(t1~t2) +
    theme(axis.text.x = element_text(angle=45,hjust = 1,vjust = 1)) +
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=12),
          strip.text = element_text(size=10)) +
    labs(x="Distance (\u03bcm)",y="Cross-correlation")
}
