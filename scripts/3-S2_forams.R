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
forams <- read.csv("datasets/S2/S2_counts.csv", sep=";")[-1] %>%
  mutate(depth=as.numeric(gsub(",", ".", gsub("\\.", "", depth))),
         depth=depth*100) #replace commas with dots for decimals
str(forams)

# Read in forams taxa names groups
nms <- read.csv("outputs/nms_old_new.csv", na.strings = "") #handle empty cells of foraminifera indicator spp

# Read in %sand
S2_sand <- read.csv("datasets/S2/S2_sand.csv", sep = ";") %>%
  mutate(sand=as.numeric(gsub(",", ".", gsub("\\.", "", sand))),
         depth=as.numeric(gsub(",", ".", gsub("\\.", "", depth))),
         depth=depth*100)
str(S2_sand)

# Read in age-depth model
ages <- read.table("Bacon_runs/S2_core_v2/S2_core_v2_104_ages.txt")
str(ages)
colnames(ages) <- ages[1,]
ages <- ages[-1,]
ages <- data.frame(apply(ages, 2, as.numeric)) #transform to numeric

# Read in XRF data
S2_geochem_data <- read.csv("datasets/S2/S2_XRF.csv") %>%
  #dplyr::rename(depth=ï..depth) %>%
  mutate(Sr_Rb=Sr/Rb) %>% #unweathered terrestrial fraction
  mutate(lnCa_Ti=log(Ca/Ti))%>% #biogenic CaCO3 vs detrital input
  mutate(Si_Zr=Si/Zr) %>%
  mutate(depth=depth*100) %>%
  left_join(ages[c("depth", "mean")], by=c("depth")) 
str(S2_geochem_data)

vec <- c(48,49,53,54,90,92,95,97,100,102,105,107)
ages_corr <- c(1422,1581,2189,2340,8315,8601,9026,9313,9738,9939,10236,10433)

S2_geochem_data$mean[vec] <- ages_corr
names(S2_geochem_data$mean)

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
  gather(key = taxa, value = count, -depth, -sample_id) %>%
  mutate(taxa= plyr::mapvalues(taxa, from=nms$old, to=nms$new)) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = nms[,1], to = nms$wall_type)) %>%
  mutate(indspp = plyr::mapvalues(taxa, from = nms[,1], to = nms$habitat)) %>%
  mutate(assemblage=recode(assemblage,
                           '1'="agglutinated",
                           '2'="hyaline",
                           '3'="porcellanous")) %>%
  group_by(depth, taxa, assemblage, indspp) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  #filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(total_sample = sum(count)) %>% 
  full_join(ages, by=c("depth")) %>%
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
  arrange(taxa) %>%
  mutate(n_valves = ifelse(total_sample < 100, "low", "high")) %>%
  filter(!is.na(indspp))#here label samples if nvalves >100 or higher
  #drop_na()

#make it wide
core_counts_wide_forams <- core_counts_common %>%
  select(depth, mean, taxa, assemblage, indspp, relative_abundance_percent) %>%
  rename(age_calyr=mean) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(depth) #sort by increasing time

is.na(core_counts_wide_forams$age_calyr)
is.na(core_counts_wide_forams$indspp)

# here assign manually age to depths that join did not work for some reason
core_counts_wide_forams[c(22,23),2] <- 1581

core_counts_wide_forams[c(57,58,59,60,61,62),2] <- 1581

core_counts_wide_forams[c(31,32),2] <- 1808
core_counts_wide_forams[c(36,37),2] <- 1961
core_counts_wide_forams[c(218,219,220),2] <- 8315
core_counts_wide_forams[c(224),2] <- 8601
core_counts_wide_forams[c(227),2] <- 8790
core_counts_wide_forams[c(229),2] <- 9026
core_counts_wide_forams[c(230,231),2] <- 9073
core_counts_wide_forams[c(243,244),2] <- 9738
core_counts_wide_forams[c(250,251),2] <- 9939
core_counts_wide_forams[c(253,254),2] <- 10073
core_counts_wide_forams[c(256,257),2] <- 10236

#do coniss to add statistically significant stratigraphic zones
core_counts_wide_forams[is.na(core_counts_wide_forams)] <- 0

core_counts_wide_forams_unique <- core_counts_wide_forams %>%
  select(-assemblage) %>%
  distinct(depth, .keep_all=TRUE)

core_counts_wide_forams_unique[is.na(core_counts_wide_forams_unique)] <- 0

foramsHel <- decostand(core_counts_wide_forams_unique[,3:ncol(core_counts_wide_forams_unique)], 
                       method="hellinger")
diss <- vegdist(foramsHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=3) #k=groups
locate <- cumsum(rle(zones)$lengths)+1
zones <- core_counts_wide_forams_unique[locate, ][,2] #by age #ULL AQUI pq tinc la matriu de dades doblada per assemblage
zones <- zones$age_calyr

## Plot S2 forams record 
#theme_set(theme_bw(12))
theme_set(theme_paleo())

stratiplot <- ggplot(core_counts_common, aes(x = relative_abundance_percent, y = mean, colour=indspp)) +
  geom_col_segsh(size=1.3) +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa), rotate_facet_labels = 70) +
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "Age (cal yrs BP)", colour="Assemblage") +
  ggtitle("S2 core") +
  theme (legend.position = "bottom") +
  theme(axis.text.x=element_text(size = 8)) +
  #geom_hline(yintercept = zones, col = "black", lty = 2, alpha = 0.9)
stratiplot

## XRF
#do coniss on XRF to add statistically significant stratigraphic zones
S2_geochem_data <- S2_geochem_data %>%
  rename(age_calyr=mean) 

S2_geochem_data[is.na(S2_geochem_data)] <- 0

XRFHel <- decostand(S2_geochem_data[,2:22], method="hellinger")
diss <- vegdist(XRFHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=3) #k=groups
locate <- cumsum(rle(zones)$lengths)+1
zones_XRF <- S2_geochem_data[locate, ][,ncol(S2_geochem_data)]
#zones[1] <- NA
#zones <- zones$depth

## Plot XRF
# Prepare long-format
S2_geochem_long <- gather(data=S2_geochem_data, key = param, value = value, -depth, -age_calyr) %>% #depth is column 1
  filter(!age_calyr==0)
  
S2_plot_geochem <- S2_geochem_long %>%
  filter(param %in% c("Ca", "Si", "S", "Fe", "Zr")) %>%
  ggplot(aes(x = value, y = age_calyr)) +
  geom_lineh() +
  geom_hline(yintercept = zones_XRF, col = "black", lty = 2, alpha = 0.9) +
  geom_point() +
  scale_y_reverse() +
  facet_geochem_gridh(vars(param)) +
  labs(x = NULL, y = "Age (cal yr BP)") +
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
  dplyr::select(-c(1)) %>% 
  gather(key = param, value = value, -depth) %>%
  left_join(ages[c("depth", "mean")], by="depth")

sand_depth <- c(66,70,72,136,139,141,143,144,149,152,154,156)
ages_corr <- c(1581,1808,1961,8315,8601,8790,9026,9073,9738,9939,10073,10236)

S2_sand_long$mean[sand_depth] <- ages_corr

S2_sand_nonmissing <- S2_sand_long %>%
  arrange(depth)

S2_sand_plt <- ggplot(S2_sand_long, aes(x = mean, y = value)) +
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

write.csv(S2_sand_long, "outputs/S2_sand_age.csv")

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

#ggsave("outputs/S2_multiproxy_new_18092023.png", plt, height = 6, width = 10)

## Composite plot (stratiplots + PCA)
# library(cowplot)
# 
# composite <- ggdraw() +
#   draw_plot(plot) +
#   draw_plot(PCA_XRF, x = 0.8, y = 0.7, width = .3, height = .3)
# composite
# 
# ggsave("outputs/composite.png", composite, height = 6, width = 10)
