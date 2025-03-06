library(tidyverse)
library(magrittr)

filename_manage <- \(x){
  strsplit(x, "/") %>% `[[`(1) %>% tail(1) %>% gsub("_cell_seg_data.txt", "", .)
}

pdac_raw <- Sys.glob("examples/data/2024Phenotypes/PDAC_ModPhenotypes/*.txt") %>%
  lapply(\(f) read.delim(f,sep=",") %>% mutate(Spot=filename_manage(f))) %>%
  bind_rows()

pdac_raw %>%
  count(Spot)

pdac_raw %>%
  filter(grepl("503", Spot)) %>%
ggplot(aes(Cell.X.Position, Cell.Y.Position, color=CellOfInterest,shape = CellOfInterest)) +
  geom_point(size=2) +
  facet_wrap(~ Spot) +
  theme_void() +
  scale_color_manual(values = as.vector(pals::glasbey(n=6)),drop=FALSE)

pdac_raw %>% duplicated() %>% sum()

ipmn_raw <- Sys.glob("examples/data/2024Phenotypes/IPMN_ModPhenotypes/*.txt") %>%
  lapply(\(f) read.delim(f,sep=",") %>% mutate(Spot=filename_manage(f))) %>%
  bind_rows()

ipmn_raw %>%
  count(CellOfInterest)

ipmn_raw %>%
  mutate(CellOfInterest = factor(CellOfInterest)) %>%
  filter(grepl("784", Spot)) %>%
  ggplot(aes(Cell.X.Position, Cell.Y.Position, color=CellOfInterest,shape = CellOfInterest)) +
  geom_point(size=2) +
  facet_wrap(~ Spot) +
  theme_void() +
  scale_color_manual(values = as.vector(pals::glasbey(n=6)),drop=FALSE)

ipmn_raw %>% duplicated() %>% sum()


markers <- c("Treg", "APC", "Epithelial", "HelperT", "PDL1_CD3", "PDL1_CD8", "PDL1_FoxP3", "CD4", "CTLs")
for(mark in markers){
  pdac_raw %<>% mutate(!!mark := ifelse(pdac_raw[,mark] == "neg", 0, 1))
}

pdac_signs <- pdac_raw %>% unite("Signature", markers, sep = "")

pdac_processed <- pdac_raw %>%
  left_join(pdac_signs) %>% # will complain! why
  mutate(NoMark = ifelse(Signature == "000000000", 1, 0)) %>%
  dplyr::select(-Signature)

df_plot <- pdac_processed %>%
  filter(grepl("513", Spot)) %>%
  pivot_longer(cols = -c(contains("Cell"), "Spot"), names_to = "Marker") %>%
  dplyr::filter(value == 1)

ggplot(df_plot, aes(Cell.X.Position, Cell.Y.Position, color=Marker)) +
  geom_point(size=.5) +
  facet_wrap(~ Spot) +
  theme_void()


# Count occurrences of join keys in `pdac_raw`
pdac_raw_dupes <- pdac_raw %>%
  count(Cell.X.Position, Cell.Y.Position, Spot) %>%
  filter(n > 1)

# Count occurrences of join keys in `pdac_signs`
pdac_signs_dupes <- pdac_signs %>%
  count(Cell.X.Position, Cell.Y.Position, Spot) %>%
  filter(n > 1)

# Identify many-to-many keys (appear >1 in both tables)
many_to_many_keys <- inner_join(pdac_raw_dupes, pdac_signs_dupes,
                                by = c("Cell.X.Position", "Cell.Y.Position", "Spot"))

# Print results
print("Duplicate keys in pdac_raw:")
print(pdac_raw_dupes)

print("Duplicate keys in pdac_signs:")
print(pdac_signs_dupes)

print("Many-to-Many keys (appear in both with duplicates):")
print(many_to_many_keys)

pdac_raw %>%
  inner_join(pdac_raw_dupes)
