library(posterior)

sim_idx <- 6
nx <- ny <- 30
n_thin <- 5
out_full <- readRDS(paste0("/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingk_",sim_idx,".rds"))
vhat <- out_full[[1]]$v_mcmc
vhat <- lapply(vhat,\(v) v[1:(nx*ny),])
arr <- simplify2array(vhat)
arr <- aperm(arr,perm=c(3,1,2))
vhat <- posterior::rvar(arr)

beta <- out_full[[1]]$beta_mcmc
arr <- beta
arr <- arr[,,seq(1,dim(arr)[3],by=n_thin),drop=FALSE]
arr <- aperm(arr,perm=c(3,1,2))
beta <- posterior::rvar(arr)
out_exp <- lapply(out_full,\(o) list(v_mcmc=E(vhat),
                                     beta_mcmc=E(beta)))

saveRDS(out_exp,paste0("/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingk_",sim_idx,"_exp.rds"))

# then paste the following to get the exp back
# cpdn /nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingk_6_exp.rds examples/data
