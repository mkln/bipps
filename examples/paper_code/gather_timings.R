sim_idx <- 1:30

l_k <- unlist(lapply(sim_idx,\(i) readRDS(paste0("timings_varyingk_",i,".rds"))))
l_px <- unlist(lapply(sim_idx,\(i) readRDS(paste0("timings_varyingpx_",i,".rds"))))

df <- data.frame(sim_idx=c(sim_idx,sim_idx),timing=c(l_k,l_px),type=c(rep("k",30),rep("px",30)))

saveRDS(df,"timings.rds")
