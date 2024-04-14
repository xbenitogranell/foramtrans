#Clear workspace
rm(list=ls(all=TRUE))

##loading libraries for functions used
library(analogue) 
library(rioja) 
library(plyr) #allow to join dataframes by common column
library(dplyr) #allow to summarise variables and manipulate multiple dataframes
library(ggplot2) #to make nice plots
library(tidyverse)
library(tidypaleo) #(https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(cluster)
library(mgcv)

## Read in forams counts
carlet <- read.csv("datasets/Carlet/carlet_forams_counts.csv", sep=",")[-1] %>%
  mutate(depth=depth*-1)
stjaume <- read.csv("datasets/Sant Jaume/stjaume_forams_counts.csv", sep=",")[-1] %>%
  mutate(depth=depth*-1)

# Read nms list to standardize
nms <- read.csv("outputs/nms_old_new.csv")

# Read in age-depth models
ages_carlet <- read.table("Bacon_runs/carlet/carlet_156_ages.txt")
str(ages_carlet)
colnames(ages_carlet) <- ages_carlet[1,]
ages_carlet <- ages_carlet[-1,]
ages_carlet <- data.frame(apply(ages_carlet, 2, as.numeric)) #transform to numeric

ages_stjaume <- read.table("Bacon_runs/stjaume_v2/stjaume_v2_146_ages.txt")
str(ages_stjaume)
colnames(ages_stjaume) <- ages_stjaume[1,]
ages_stjaume <- ages_stjaume[-1,]
ages_stjaume <- data.frame(apply(ages_stjaume, 2, as.numeric)) #transform to numeric

# Tweak the dataset
new <- stjaume %>% 
  gather(key = taxa, value = count, -depth) %>%
  mutate(taxa=plyr::mapvalues(taxa, from=nms$old, to=nms$new)) %>%
  mutate(assemblage = plyr::mapvalues(taxa, from = nms[,1], to = nms$wall_type)) %>%
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
  mutate(total_sample = sum(count)) %>% 
  left_join(ages_stjaume[c("depth", "min", "max", "mean")], by=c("depth")) %>%
  ungroup()

levels(factor(new$taxa))
levels(factor(new$assemblage))

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
  mutate(n_valves = ifelse(total_sample < 100, "low", "high")) #here label samples if nvalves >100 or higher
#drop_na()

#make it wide
core_counts_wide_forams <- core_counts_common %>%
  select(depth, mean, min,max, taxa, assemblage, count) %>%
  rename(age_calyr=mean) %>%
  rename(upper_age=min) %>%
  rename(lower_age=max) %>%
  spread(key = taxa, value = count) %>%
  arrange(depth) #sort by increasing time

is.na(core_counts_wide_forams$age_calyr)
median(diff(core_counts_wide_forams$upper_age),na.rm = TRUE)

# assign df to each core
#stjaume_counts_wide <- core_counts_wide_forams
#write.csv(stjaume_counts_wide, "datasets/Sant Jaume/stjaume_counts_wide.csv")

stjaume_counts_wide <- core_counts_wide_forams
#write.csv(carlet_counts_wide, "datasets/Carlet/carlet_counts_wide.csv")

##
# Read in %sand
stjaume_sand <- read.csv("datasets/Sant Jaume/stjaume_sand.csv", sep = ",") %>%
  mutate(depth=depth*-1) %>%
  left_join(ages_stjaume[c("depth", "mean")], by=c("depth")) 

str(stjaume_sand)

#do coniss to add statistically significant stratigraphic zones
core_counts_wide_forams[is.na(core_counts_wide_forams)] <- 0

core_counts_wide_forams_unique <- core_counts_wide_forams %>%
  select(-assemblage) %>%
  distinct(depth, .keep_all=TRUE)

core_counts_wide_forams_unique[is.na(core_counts_wide_forams_unique)] <- 0

foramsHel <- decostand(core_counts_wide_forams_unique[,4:ncol(core_counts_wide_forams_unique)], 
                       method="hellinger")
diss <- vegdist(foramsHel, method="bray")
clust <- chclust(diss, method="coniss")
bstick(clust)

zones <- cutree(clust, k=4) #k=groups
locate <- cumsum(rle(zones)$lengths)+1
zones <- core_counts_wide_forams_unique[locate, ][,2] #by age #ULL AQUI pq tinc la matriu de dades doblada per assemblage
zones <- zones$age_calyr

## Plot forams record 
#theme_set(theme_bw(12))
theme_set(theme_paleo())

stratiplot <- ggplot(core_counts_common, aes(x = relative_abundance_percent, y = mean, colour=assemblage)) +
  geom_col_segsh(size=1.3) +
  scale_y_reverse() +
  facet_abundanceh(vars(taxa), rotate_facet_labels = 70) +
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2) +
  labs(x = "Relative abundance (%)", y = "Age (cal yrs BP)", colour="Assemblage") +
  ggtitle("St Jaume core") +
  theme (legend.position = "bottom") +
  theme(axis.text.x=element_text(size = 8)) +
  geom_hline(yintercept = zones, col = "black", lty = 2, alpha = 0.9)
stratiplot

## Plot Proportion of sand
stjaume_sand_plt <- ggplot(stjaume_sand, aes(x = mean, y = sand)) +
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
stjaume_sand_plt

## Plot foram record + % sand
library(patchwork)

plt <- wrap_plots(
  stratiplot + 
    theme(strip.background = element_blank(), strip.text.y = element_blank()),
  stjaume_sand_plt +
    theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  # PCA_components_plt +
  #   theme(axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()),
  nrow = 1,
  widths = c(3,1,0.5)
)
plt

ggsave("outputs/StJaume.png", plt, height = 6, width = 10)
