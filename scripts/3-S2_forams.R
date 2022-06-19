## Exploratory analyses of S2 Foraminiferal record (Ebro Delta, Spain)
# author: Xavier Benito (xavier.benito.granell@gmail.com)'
# date: 'Date: 18/06/2022'

##loading libraries for functions used
library(analogue) #to join diatom datasets on their common spp
library(rioja) #to merge diatom datasets on their common spp
library(plyr) #allow to join dataframes by common column
library(dplyr) #allow to summarise variables and manipulate multiple dataframes
library(ggplot2) #to make nice plots
library(tidyverse)
library(tidypaleo) #(https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(cluster)

## Read in forams counts
forams <- read.csv("datasets/S2/S2_counts.csv") 
str(forams)

# Read in forams taxa names groups
changes <- read.csv("datasets/S2/nms_taxa_groups.csv")

# Read in XRF data
S2_geochem_data <- read.csv("datasets/S2/S2_XRF.csv")
str(S2_geochem_data)

## Calculate relative abundance
# sample_info <- forams[, names(forams) %in% c("section", "sample_id", "depth")] 
# forams <- forams[, !names(forams) %in% c("section", "sample_id", "depth")]
# forams[is.na(forams)] <- 0

# Transform to relative abundance
# total <- apply(forams, 1, sum)
# forams <- forams/total*100
# 
# ##Remove rare species
# abund <- apply(forams, 2, max)
# n.occur <- apply(forams>0, 2, sum)
# forams <- forams[, n.occur>1 & abund>2] #more than 2% of RA and present in >1 sample


#this is to transform to tidy format, calculate % and subset more common species
new <- forams %>% 
  dplyr::select(-section) %>%
  gather(key = taxa, value = count, -depth, -sample_id) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = changes$taxa, to = changes$wall_type)) %>%
  mutate(assemblage=recode(assemblage,
                           '1'="agglutinated",
                           '2'="hyaline",
                           '3'="porcellanous")) %>%
  group_by(depth, taxa, assemblage) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  #filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  ungroup()


# filter more abundant taxa; format need to be on long-wide format-->no spreaded 
core_common_taxa <- new %>%
  group_by(taxa) %>%
  summarise(max_rel_abund = max(relative_abundance_percent)) %>%
  filter(max_rel_abund >= 4) %>%
  arrange(max_rel_abund) %>%
  pull(taxa)

# select from initial table
core_counts_common <- new %>%
  filter(taxa %in% core_common_taxa) %>%
  mutate(taxa = factor(taxa, levels = core_common_taxa)) %>%
  arrange(taxa)

#make it wide
core_counts_wide_forams <- core_counts_common %>%
  select(depth, taxa, assemblage, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(depth) #sort by increasing time

#do coniss to add statistically significant stratigraphic zones
core_counts_wide_forams[is.na(core_counts_wide_forams)] <- 0

foramsHel <- decostand(core_counts_wide_forams[,3:ncol(core_counts_wide_forams)], method="hellinger")
diss <- vegdist(foramsHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=3)
locate <- cumsum(rle(zones)$lengths)+1
zones <- core_counts_wide_forams[locate, ][,1]
zones <- zones$depth

## Plot S2 forams record 
theme_set(theme_bw(12))
#theme_set(theme_paleo())

stratiplot <- ggplot(core_counts_common, aes(x = relative_abundance_percent, y = depth, colour=assemblage)) +
  geom_col_segsh(size=1.3) +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa), rotate_facet_labels = 70) +
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "core depth (m)", colour="Assemblage") +
  ggtitle("S2 core") +
  theme (legend.position = "bottom") +
  geom_hline(yintercept = zones, col = "blue", lty = 1, alpha = 0.7)
stratiplot

## Plot XRF
# Prepare long-formt
S2_geochem_long <- gather(data=S2_geochem_data, key = param, value = value, -depth)

S2_plot_geochem <- S2_geochem_long %>%
  #filter(param %in% c("Fe", "K", "Si")) %>%
  ggplot(aes(x = value, y = depth)) +
  geom_lineh() +
  geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)")
S2_plot_geochem

#ggsave("outputs/forams_stratplot.png", stratiplot, height = 6, width = 10)

# Adding dendrograms in geochem data
coniss <- S2_geochem_long %>%
  nested_data(qualifiers = c(depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

S2_plot_geochem +
  layer_dendrogram(coniss, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss, aes(y = depth))

# Combine species abundance data and XRF plots
library(patchwork)
wrap_plots(
  stratiplot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  S2_plot_geochem +
    #layer_dendrogram(coniss, component = "CONISS", aes(y = depth)) +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  nrow = 1,
  widths = c(1, 1)
)
