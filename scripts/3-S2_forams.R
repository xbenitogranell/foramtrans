## Exploratory analyses of S2 Foraminiferal record (Ebro Delta, Spain)
# author: Xavier Benito (xavier.benito.granell@gmail.com)'
# date: 'Date: 18/06/2022'

#Clear workspace
rm(list=ls(all=TRUE))

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

# Read in %sand
S2_sand <- read.csv("datasets/S2/S2_sand.csv", sep = ";")
str(S2_sand)

# replace commas to points
S2_sand$depth <- as.numeric(gsub(",", ".", gsub("\\.", "", S2_sand$depth))) #replace commas with dots for decimals
S2_sand$sand <- as.numeric(gsub(",", ".", gsub("\\.", "", S2_sand$sand))) #replace commas with dots for decimals

# Read in XRF data
S2_geochem_data <- read.csv("datasets/S2/S2_XRF.csv") %>%
  dplyr::rename(depth=ï..depth) %>%
  mutate(Sr_Rb=Sr/Rb) %>% #unweathered terrestrial fraction
  mutate(lnCa_Ti=log(Ca/Ti))%>% #biogenic CaCO3 vs detrital input
  mutate(Si_Zr=Si/Zr) 
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
  dplyr::select(-1) %>%
  gather(key = taxa, value = count, -depth, -sample_id) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = changes[,1], to = changes$wall_type)) %>%
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
  filter(max_rel_abund >= 10) %>%
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
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "core depth (m)", colour="Assemblage") +
  ggtitle("S2 core") +
  theme (legend.position = "bottom") +
  theme(axis.text.x=element_text(size = 8)) +
  geom_hline(yintercept = zones, col = "black", lty = 2, alpha = 0.9)
stratiplot


## XRF
#do coniss on XRF to add statistically significant stratigraphic zones
S2_geochem_data[is.na(S2_geochem_data)] <- 0

XRFHel <- decostand(S2_geochem_data[,2:ncol(S2_geochem_data)], method="hellinger")
diss <- vegdist(XRFHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=3) #k=groups
locate <- cumsum(rle(zones)$lengths)+1
zones <- S2_geochem_data[locate, ][,1]
#zones <- zones$depth

## Plot XRF
# Prepare long-format
S2_geochem_long <- gather(data=S2_geochem_data, key = param, value = value, -depth) #depth is column 1

S2_plot_geochem <- S2_geochem_long %>%
  filter(param %in% c("Ca", "Si", "S", "Fe", "Zr")) %>%
  ggplot(aes(x = value, y = depth)) +
  geom_lineh() +
  geom_hline(yintercept = zones, col = "black", lty = 2, alpha = 0.9) +
  #geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Depth (cm)") +
  theme(axis.text.x=element_text(size = 8))
S2_plot_geochem

#ggsave("outputs/forams_stratplot.png", stratiplot, height = 6, width = 10)

##Adding dendrograms in geochem data
coniss <- S2_geochem_long %>%
  nested_data(qualifiers = c(depth), key = param, value = value, trans = scale) %>%
  nested_chclust_coniss()

S2_plot_geochem +
  layer_dendrogram(coniss, aes(y = depth), param = "CONISS") +
  layer_zone_boundaries(coniss, aes(y = depth))

## Plot Proportion of sand
S2_sand_long <- S2_sand %>%
  dplyr::select(-c("ï..section","sample_id")) %>% 
  gather(key = param, value = value, -depth)

S2_sand_plt <- ggplot(S2_sand_long, aes(x = depth, y = value)) +
  # geom_areah() +
  geom_line() +
  # geom_point() +
  scale_y_reverse() +
  scale_x_reverse() +
  #facet_grid(~param) +
  coord_flip() +
  theme(axis.text.x=element_text(size = 8)) +
  #geom_smooth() +
  labs(x = NULL, y = "sand (%)")
  #ggtitle("S2 % sand")
S2_sand_plt

## Plot XRF elements + % sand
library(patchwork)

plt <- wrap_plots(
  S2_plot_geochem +
    #layer_dendrogram(coniss, component = "CONISS", aes(y = depth)) +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  S2_sand_plt +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  nrow = 1,
  widths = c(3,0.5,1)
)
plt


# Combine species abundance data and XRF plots
plt <- wrap_plots(
  stratiplot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  S2_plot_geochem +
    #layer_dendrogram(coniss, component = "CONISS", aes(y = depth)) +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    labs(y = NULL),
  S2_sand_plt +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  # PCA_components_plt +
  #   theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  nrow = 1,
  widths = c(3,1,0.5)
)
plt


#ggsave("outputs/S2_multiproxy.png", plt, height = 6, width = 10)

## Composite plot (stratiplots + PCA)
# library(cowplot)
# 
# composite <- ggdraw() +
#   draw_plot(plot) +
#   draw_plot(PCA_XRF, x = 0.8, y = 0.7, width = .3, height = .3)
# composite
# 
# ggsave("outputs/composite.png", composite, height = 6, width = 10)
