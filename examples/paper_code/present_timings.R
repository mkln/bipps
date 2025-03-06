sim_timings <- readRDS("examples/data/timings.rds")

real_timings <- readRDS("examples/data/timings_apps.rds")

szs <- c(12,24,48)
# mcmc and model settings
grid <- expand.grid(sz=szs,sim=1:10) %>%
  mutate(sim_idx = 1:nrow(.)) %>%
  mutate(sz=as.factor(sz))

sim_timings %>%
  filter(type == "px") %>%
  left_join(grid) %>%
  ggplot(aes(sz,timing)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  labs(x="Number of pixels",y="Model fitting time (s)")
