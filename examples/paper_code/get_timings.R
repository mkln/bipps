
args <- commandArgs(trailingOnly = TRUE)
sim_idx <- as.integer(args[1])

out_full <- readRDS(paste0("/nfs/turbo/umms-ukarvind/joelne/bipps_simulations/out_simset_varyingk_",sim_idx,".rds"))

timing <- unlist(lapply(out_full,\(o) o$mcmc_time))

saveRDS(timing,paste0("timings_varyingk_",sim_idx,".rds"))
