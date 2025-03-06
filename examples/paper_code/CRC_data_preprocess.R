# from original CRC dataset

df_raw <- read_csv("data/CRC_master.csv") %>%
  mutate(type = as.factor(type)) %>%
  filter(!(type %in% c("dirt","undefined"))) %>%
  mutate(type = fct_collapse(type,
                             stroma = c("stroma","nerves","lymphatics"),
                             "CD163+ macros" = c("CD68+CD163+ macrophages","CD163+ macrophages"),
                             "CD68+ macros" = c("CD11b+ monocytes","CD11b+CD68+ macrophages","CD68+ macrophages",
                                                "CD68+ macrophages GzmB+","CD11c+ DCs"),
                             "generic immune" = c("tumor cells / immune cells","immune cells","immune cells / vasculature"),
                             "memory CD4+ T"="CD4+ T cells CD45RO+",
                             "CD4+ T cells" = c("CD4+ T cells","CD4+ T cells GATA3+","CD3+ T cells"))) %>%
  droplevels()
