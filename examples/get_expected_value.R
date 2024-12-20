library(posterior)

sim_idx <- 6
out_full <- readRDS(paste0("/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset1_group_diff_",sim_idx,".rds"))
out_exp <- lapply(out_full,\(o) list(v_mcmc=E(o$v_mcmc),
                                     beta_mcmc=E(o$beta_mcmc)))

saveRDS(out_exp,paste0("/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset1_group_diff_",sim_idx,"_exp.rds"))
