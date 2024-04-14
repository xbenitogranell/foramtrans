# Timetrack analyses comparing forams training set and cores

#clear workspace
rm(list=ls(all=TRUE))
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

##loading libraries for functions used
library(analogue) #to perform timetrack analysis
library(tidyverse) #to manipulate dataframes
library(rioja) #to perform constrained hierarchical clustering
library(mgcv)
library(cluster)
library(ggplot2)
library(viridis)

# Load functions needed--timetack
source("scripts/functions.R")

## Read in forams counts
carlet <- read.csv("datasets/Carlet/carlet_forams_counts.csv", sep=",")[-1] %>%
  mutate(depth=depth*-1)
stjaume <- read.csv("datasets/Sant Jaume/stjaume_forams_counts.csv", sep=",")[-1] %>%
  mutate(depth=depth*-1)

# Read nms list to standardize taxa names
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
new <- carlet %>% 
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
  filter(max_rel_abund >= 2) %>%
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
  select(depth, mean, min,max, taxa, count) %>%
  rename(age_calyr=mean) %>%
  rename(upper_age=min) %>%
  rename(lower_age=max) %>%
  spread(key = taxa, value = count) %>%
  arrange(depth) #sort by increasing time

is.na(core_counts_wide_forams$age_calyr)
median(diff(core_counts_wide_forams$upper_age),na.rm = TRUE)

# assign df to each core
stjaume_counts_wide <- core_counts_wide_forams
carlet_counts_wide <- core_counts_wide_forams

## Read S2 foraminifera records
forams <- read.csv("datasets/S2/S2_counts.csv", sep=";")[-1] %>%
  mutate(depth=as.numeric(gsub(",", ".", gsub("\\.", "", depth))),
         depth=depth*100) #replace commas with dots for decimals
str(forams)

# Read in S2 age-depth model
ages <- read.table("Bacon_runs/S2_core_v2/S2_core_v2_104_ages.txt")
str(ages)
colnames(ages) <- ages[1,]
ages <- ages[-1,]
ages <- data.frame(apply(ages, 2, as.numeric)) #transform to numeric

#this is to transform to tidy format, calculate % and subset more common species
new <- forams %>% 
  gather(key = taxa, value = count, -depth, -sample_id) %>%
  mutate(taxa= plyr::mapvalues(taxa, from=nms$old, to=nms$new)) %>%
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
  left_join(ages[c("depth", "min", "max", "mean")], by=c("depth")) %>%
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
  select(depth, mean, min,max, taxa, count) %>%
  rename(age_calyr=mean) %>%
  rename(upper_age=min) %>%
  rename(lower_age=max) %>%
  spread(key = taxa, value = count) %>%
  arrange(depth) #sort by increasing time

is.na(core_counts_wide_forams$age_calyr)
median(diff(core_counts_wide_forams$upper_age),na.rm = TRUE)

# here assign manually age to depths that join did not work for some reason
core_counts_wide_forams[12,2] <- c(1581)
core_counts_wide_forams[12,3] <- c(1351)
core_counts_wide_forams[12,4] <- c(1785)

core_counts_wide_forams[16,2] <- c(1808)
core_counts_wide_forams[16,3] <- c(1589)
core_counts_wide_forams[16,4] <- c(2051)

core_counts_wide_forams[18,2] <- c(1961)
core_counts_wide_forams[18,3] <- c(1721)
core_counts_wide_forams[18,4] <- c(2240)

core_counts_wide_forams[82,2] <- c(8315)
core_counts_wide_forams[82,3] <- c(8089)
core_counts_wide_forams[82,4] <- c(8601)

core_counts_wide_forams[85,2] <- c(8601)
core_counts_wide_forams[85,3] <- c(8295)
core_counts_wide_forams[85,4] <- c(8970)

core_counts_wide_forams[87,2] <- c(8790)
core_counts_wide_forams[87,3] <- c(8445)
core_counts_wide_forams[87,4] <- c(9190)

core_counts_wide_forams[89,2] <- c(9026)
core_counts_wide_forams[89,3] <- c(8664)
core_counts_wide_forams[89,4] <- c(9431)

core_counts_wide_forams[90,2] <- c(9073)
core_counts_wide_forams[90,3] <- c(8700)
core_counts_wide_forams[90,4] <- c(9500)

core_counts_wide_forams[95,2] <- c(9738)
core_counts_wide_forams[95,3] <- c(9364)
core_counts_wide_forams[95,4] <- c(10202)

core_counts_wide_forams[98,2] <- c(9939)
core_counts_wide_forams[98,3] <- c(9586)
core_counts_wide_forams[98,4] <- c(10354)

core_counts_wide_forams[100,2] <- c(10073)
core_counts_wide_forams[100,3] <- c(9737)
core_counts_wide_forams[100,4] <- c(10483)

core_counts_wide_forams[102,2] <- c(10236)
core_counts_wide_forams[102,3] <- c(9935)
core_counts_wide_forams[102,4] <- c(10635)

S2_counts_wide <- core_counts_wide_forams

# S2 <- read.csv("datasets/S2/S2_counts.csv", sep=";")[-1]
# carlet <- read.csv("datasets/Carlet/carlet_forams_counts.csv", sep=",")[-1]
# stjaume <- read.csv("datasets/Sant Jaume/stjaume_forams_counts.csv", sep=",")[-1]

# Read in foraminifera training set
code_spp <- read.csv("datasets/training set/code_spp.csv")
head(code_spp)
code_spp$spp <- gsub(" ",".",code_spp$spp)

forams_surface <- read.csv("datasets/training set/forams_surface.csv", sep=",")
forams_surface_spp <- forams_surface[-c(1:2),]
nms_modern <- forams_surface_spp[1,]
colnames(forams_surface_spp) <- nms_modern[1,]
forams_surface_spp <- forams_surface_spp[-1,]
forams_surface_spp <- forams_surface_spp[,-c(156,157)]
names(forams_surface_spp)
str(forams_surface_spp)

# Transform numeric dataframe
mydata_num <- data.frame(apply(forams_surface_spp[,5:ncol(forams_surface_spp)], 2, as.numeric)) 

mydata_num[is.na(mydata_num)] <- 0
str(mydata_num)
colSums(mydata_num)
rowSums(mydata_num)

forams_surface <- forams_surface[-c(1:3),-c(156,157)]
colnames(forams_surface) <- nms_modern[-c(156,157)]

# join sites, habitat and month columns
forams_surface_df <- mydata_num %>%
  bind_cols(forams_surface[,c(1:4)]) %>%
  gather(spp, abund, -Sample, -Site,-Habitat,-Period) %>%
  mutate(spp=plyr::mapvalues(spp, from=code_spp$code, to=code_spp$spp)) %>%
  mutate(spp=plyr::mapvalues(spp, from=nms$old, to=nms$new)) %>%
  group_by(Sample, Site, Habitat, Period, spp) %>% #there are >2 samples for the same site, habitat and period
  summarise(abund = sum(abund)) %>% #
  #filter(!abund == 0) %>% #this is to remove empty samples (rows)
  spread(key = spp, value = abund) %>%
  as.data.frame()

# Make a common df from modern and fossil foraminiferal counts
df <- analogue::join(S2_counts_wide,
                     carlet_counts_wide,
                     stjaume_counts_wide,
                     forams_surface_df,
                     verbose = TRUE)

# name the list
names(df) <- c("S2", "Carlet", "StJaume", "trainingset")
spp <- colnames(df[[1]])
#write.csv(spp, "outputs/spp_all_nms.csv")

#check NA in the list
listnans <- lapply(df, function(x) sum(is.na(x)))

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "age_calyr", "upper_age", "lower_age", "Sample", "Site", "Habitat", "Period"))] # drop year & depths vars
  #core <- core[, -which(names(core) %in% c("depth", "age_calyr", "lower_age", "Sample", "Site", "Habitat", "Period"))] # drop year & depths vars
  core <- core[, colSums(core) > 0] #select only present species
  #core <- core[core$upper_age > 0,] #subset samples with age>0
  #rownames(core) <- core$upper_age
  #core <- core[, -which(names(core) %in% c("upper_age"))] # drop remaining column not wanted
  # core <- tran(core, "hellinger")
  return(core)
}

cores <- lapply(seq_along(df), remove, cores=df)
names(cores) <- names(df)

training <- cores$trainingset
cores$trainingset <- NULL

#saveRDS(cores, "outputs/forams_cores_list.rds")

#create list of metadata
mymetadata <- function(i, cores, ...) {
  core <- cores[[i]]
  rownms <- df[[i]][["sample_id"]]
  core <- core[, which(names(core) %in% c("age_calyr", "depth", "upper_age","lower_age"))] 
  rownames(core) <- rownms
  return(core)
}

metadata_cores <- lapply(seq_along(df), mymetadata, cores=df)
names(metadata_cores) <- names(df)
metadata_cores[[4]] <- NULL

#saveRDS(metadata_cores, 'outputs/agedepthForams.rds')

#Calculate SQchord dissimilarity from successive samples
#drop trainingset from the core list; first save the df 
df$trainingset <- NULL

doSQchord <- function(i, cores,..) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "age_calyr", "upper_age", "lower_age", "Sample", "Site", "Habitat", "Period"))] # drop year & depths vars
  core <- core[, colSums(core) > 0] #select only present species
  D <- analogue::distance(core/100, method="SQchord")
  # create an indicator for all diagonals in the matrix
  d <- row(D) - col(D)
  # use split to group on these values
  diagonal<-split(D, d)
  SQchord<-unlist(diagonal["1"]) # select relevant one (diag = 1)
  upper_age <- df[[i]]$upper_age
  lower_age <- df[[i]]$lower_age
  cbind.data.frame(SQchord, upper_age[-1], lower_age[-1]) 
}

#Wrap-up function over the cores list and name the cores
coresSQchord <- lapply(seq_along(df), doSQchord, cores=df)
names(coresSQchord) <- names(cores)

#Extract core dataframes from the list
ts_SQchord <- plyr::ldply(coresSQchord, data.frame)
colnames(ts_SQchord) <- c("core", "SQchord", "upper_age", "lower_age")

ggplot(ts_SQchord, aes(x=upper_age, y=SQchord)) +
  geom_point() + geom_line() +
  facet_grid(core~., scales = "free") +
  theme_bw()

#write.csv(ts_SQchord, "outputs/ts_forams_SQchord.csv")

## Read in modern environnmental dataset
env_data <- read.csv("datasets/training set/env_var_surface.csv")[-1,]
str(env_data)
rownames(env_data) <- env_data$Sample

mydata_num <- data.frame(apply(env_data[,5:ncol(env_data)], 2, as.numeric)) 
rownames(mydata_num) <- rownames(env_data)
env_data_df <- data.frame(env_data[, which(names(env_data) %in% c("Sample", "Site", "Habitat", "Period"))],
                          mydata_num)
str(env_data_df)
forams_env <- merge(forams_surface_df,env_data_df, by="Sample")
env_surf <- forams_env[, -which(names(forams_env) %in% names(forams_surface_df))][-c(1:6)]

rownames(env_surf) <- forams_env$Sample

##Perform MULTIPLE time tracks
# this is function to calculate relative abundance from counts data
RA <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- tran(core, method="percent")
  return(core)
}

cores_ra <- lapply(seq_along(cores), RA, cores=cores)

#tran <- c("sqrt")
tran <- c("hellinger")
#tran <- NULL

doTimeTrack <- function(core, train, scale = TRUE, scaling = 3) {
  timetrack(train, core, method = "rda", scale = scale, transform = tran,
            scaling = scaling, na.rm=TRUE)
}

# Ull pq training is RA and cores are absolute counts
ttracks <- lapply(cores_ra, doTimeTrack, train = training,
                  scale = FALSE, scaling = 3)

names(ttracks) <- names(cores)

# plot timetracks using Gavin's function
par(mfrow=c(2,2))
par(mar=c(3,3,2,3))

#surf <- NULL
tt.lwd <- 1

## fit univariate surface - only once
surf <- ordisurf(ttracks[[1]]$ord, env_surf[,"Sand"],
                 method = "REML", scaling = ttracks[[1]]$scaling,
                 select = TRUE, plot = FALSE, knots = 20)

for(i in seq_along(ttracks)) {
  #png("ecuador_lakes_timetrack.png", width =12, height=6, units="in", res=300)
  ttplot(ttracks[[i]], surf = surf, lwd = tt.lwd)
  title(names(ttracks[i]))
  legend("topleft", col="#3465a4", lty=1, cex=1,
         legend= c("Sand %"), bty='n')
}

summary(surf)
sand.fitted <- surf$fitted.values

## Perform MULTIPLE constrained cluster analyses
doCluster <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- decostand(core, method="hellinger") #Hellinger transform relative abundance data
  diss <- vegan::vegdist(core, method="bray")
  clust <- chclust(diss, method="coniss")
  #bstick(clust)
  clust
}

coreCluster <- lapply(seq_along(cores_ra), doCluster, cores = cores_ra)
names(coreCluster) <- names(cores)

## plot Clusters
par(mfrow=c(3,3))
par(mar=c(4,2,2,2))

# Estimate broken stick model for significant zones (visual inspection)
for (i in seq_along(coreCluster)) {
  core <- coreCluster[[i]]
  n <- ncol(cores[[i]])
  k <- seq(1,n)
  sumFrac <- 1/k
  bstick <- rep(NA,n)
  for(j in 1:n) bstick[j] = 1/n*sum(sumFrac[j:n])
  
  # Compare variance explained by each split
  clustVarEx <- rev(diff(core$height)/max(core$height))
  plot(k[1:10], clustVarEx[1:10], type="o", col="black", xlab="Number of Splits", ylab="% Var Expl")
  title(names(coreCluster[i]))
  points(k[1:10], bstick[1:10], type="o", col="red")
  
}

# Estimate broken stick model for significant zones (visual inspection)
for (i in seq_along(coreCluster)) {
  bstick(coreCluster[[i]])
  title(names(coreCluster[i]))
}

#extract core clusters
nams <- names(coreCluster)
for (i in seq_along(coreCluster)) {
  assign(paste0("cCluster_", nams[i]), coreCluster[[i]])
}

#Extract significant diatom zones
# k=number of groups, so split the groups of the ccluster analysis accordingly
sigClustS2 <- cutree(coreCluster$S2, k=3)
locate <- cumsum(rle(sigClustS2)$lengths)+1
zones_S2 <- df$S2$upper_age[locate]

#write.csv(sigClustS2,"outputs/sigClustS2.csv")

sigClustCarlet <- cutree(coreCluster$Carlet, k=3)
locate <- cumsum(rle(sigClustCarlet)$lengths)+1
zones_Carlet <- df$Carlet$upper_age[locate]

## Make multittplot     
png("outputs/phase-space.png",
    res = 300,
    width = 12, height = 10, units = 'in')

par(mfrow=c(4,5))
par(mar=c(2,2,2,2))

core <- c("S2")
core_ages <- c("S2")
k=max(unique(sigClustS2)) #number of groups in ccluster

tt.lwd <- 1
env_data <- env_surf %>%
  select(Water.depth, Sand, Salinity, OM)

ordResult <- list()

for(i in 1:length(env_data)){
  tpch = 21
  bpch = 22
  pcol = "#a40000"
  pbg = "#ef2929"
  pcex = 1.6
  
  apch = 24
  acol = "orange"
  abg = "orange"
  
  multittplot(ttracks[[core]], lwd = tt.lwd, surf=NULL)
  
  coreCluster_site <- coreCluster[[core]]
  sigClust <- cutree(coreCluster_site, k=k)
  
  agesClust <- as.numeric(df[[core_ages]]$upper_age)
  depths <- as.numeric(df[[core_ages]]$depth)

  fits <- as.data.frame(cbind(fitted(ttracks[[core]], type = "passive"), 
                              traj= sigClust, depth=depths, brks = agesClust))
  
  # 
  sand_core <- S2_sand$sand[match(fits$depth, S2_sand$depth)]
  Ca_XRF <- S2_geochem_data$Ca[match(fits$depth, S2_geochem_data$depth)]
  Zr_XRF <- S2_geochem_data$Zr[match(fits$depth, S2_geochem_data$depth)]
  Fe_XRF <- S2_geochem_data$Fe[match(fits$depth, S2_geochem_data$depth)]
  
  surf <- ordisurf(ttracks[[1]]$ord, env_data[,i], method = "REML", scaling = ttracks[[1]]$scaling,
                   select = TRUE, add = TRUE, main = colnames(env_data[i]), knots = 20, 
                   npoints = nrow(env_data), col = "#3c73ad") #npoints is the number of locations at which to evaluate the fitted surface
  
  
  segments(fits$PC1[-length(fits$PC1)],
           fits$PC2[-length(fits$PC2)],
           fits$PC1[-1L],fits$PC2[-1L], col=fits$traj, type="b")
  
  locate <- cumsum(rle(fits$traj)$lengths)+1
  points(fits[locate, ], pch=24, col=acol, bg=abg)
  text(fits[locate, ], labels=fits$brks[locate], cex = 1, offset = 0.5, pos = 4)
  
  ordResult$fits <- fits
  ordResult$dev.expl[[i]] <- summary(surf)$dev.expl
  ordResult$fitted.values[[i]] <- fitted(surf)
  ordResult$grid[[i]] <- surf$grid[1:3] 
  ordResult$s.pv[[i]] <- summary(surf)[["s.pv"]]
  ordResult$ages <- agesClust
  ordResult$depth <- depths
  
  ordResult$sand_core <- sand_core
  ordResult$Zr <- Zr_XRF
  ordResult$Ca <- Ca_XRF
  ordResult$Fe <- Fe_XRF
  
  # Make prediction of core samples (PCA scrs)
  newdat <- fits[,1:2]
  
  #calibrate () calls predict.gam
  pred_space <- calibrate(surf, newdata = newdat, se.fit=TRUE) #
  
  ordResult$modspace_fit[[i]] <- as.numeric(pred_space$fit)
  ordResult$modspace_fit_error[[i]] <- as.numeric(pred_space$se.fit)
  
  # example from Paul Zander
  ordResult$RMSE[[i]] <- sqrt(mean((env_data[,i] - ordResult$modspace_fit[[i]])^2, na.action="na.exclude"))
  ordResult$split_RMSE[[i]] <- paste(round(ordResult$RMSE[[i]], digits = 2), " (", round(100 * ordResult$RMSE[[i]] / (max(env_data[,i]) - min(env_data[,i])), digits = 2), "%)")
  
  ordResult$cor[[i]] <- cor.test(sand_core, ordResult$modspace_fit[[i]], method="spearman")
  plot(sand_core,ordResult$modspace_fit[[i]], col=as.factor(ordResult$fits$traj), pch=20,
       main = paste0("Commmunity inferred-",colnames(env_data[i])," \n vs %sand"))
  legend("bottomright",legend = paste0("rho=",round(ordResult$cor[[i]]$estimate, digits=2), "\n p=", round(ordResult$cor[[i]]$p.value, digits=2)),
         inset=0.010,cex=0.9, bg="white", bty="n")
  
  ordResult$cor[[i]] <- cor.test(Zr_XRF, ordResult$modspace_fit[[i]], method="spearman")
  plot(Zr_XRF,ordResult$modspace_fit[[i]], col=as.factor(ordResult$fits$traj), pch=20,
       main = paste0("Commmunity inferred-",colnames(env_data[i])," \n vs Zr"))
  legend("bottomright",legend = paste0("rho=",round(ordResult$cor[[i]]$estimate, digits=2), "\n p=", round(ordResult$cor[[i]]$p.value, digits=2)),
         inset=0.010,cex=0.9, bg="white", bty="n")
  
  ordResult$cor[[i]] <- cor.test(Ca_XRF, ordResult$modspace_fit[[i]], method="spearman")
  plot(Ca_XRF,ordResult$modspace_fit[[i]], col=as.factor(ordResult$fits$traj), pch=20,
       main = paste0("Commmunity inferred-",colnames(env_data[i])," \n vs Ca"))
  legend("bottomright",legend = paste0("rho=",round(ordResult$cor[[i]]$estimate, digits=2), "\n p=", round(ordResult$cor[[i]]$p.value, digits=2)),
         inset=0.010,cex=0.9, bg="white", bty="n")
  
  ordResult$cor[[i]] <- cor.test(Fe_XRF, ordResult$modspace_fit[[i]], method="spearman")
  plot(Fe_XRF,ordResult$modspace_fit[[i]], col=as.factor(ordResult$fits$traj), pch=20,
       main = paste0("Commmunity inferred-",colnames(env_data[i])," \n vs Fe"))
  legend("bottomright",legend = paste0("rho=",round(ordResult$cor[[i]]$estimate, digits=2), "\n p=", round(ordResult$cor[[i]]$p.value, digits=2)),
         inset=0.010,cex=0.9, bg="white", bty="n")
  
  # 
  # abline(lm(toPlot[,1]~toPlot[,2], na.action = "na.exclude"))
  # 
  # abline(lm(sand_core~ordResult$modspace_fit[[i]], na.action="na.exclude"))
  # abline(lm(Zr_XRF~ordResult$modspace_fit[[i]], na.action="na.exclude"))
  # abline(lm(Ca_XRF~ordResult$modspace_fit[[i]], na.action="na.exclude"))
  # abline(lm(Fe_XRF~ordResult$modspace_fit[[i]], na.action="na.exclude"))

}
dev.off()

#name ordResult list  
names(ordResult$fitted.values) <- colnames(env_data)
names(ordResult$dev.expl) <- colnames(env_data)
names(ordResult$s.pv) <- colnames(env_data)
names(ordResult$grid) <- colnames(env_data)
names(ordResult$modspace_fit) <- colnames(env_data)
names(ordResult$modspace_fit_error) <- colnames(env_data)
names(ordResult$RMSE) <- colnames(env_data)

ordResult$dev.expl <- ordResult$dev.expl[order(unlist(ordResult$dev.expl),decreasing=TRUE)]
ordResult$RMSE <- ordResult$RMSE[order(unlist(ordResult$RMSE),decreasing = TRUE)]

# plot % deviance explained
for (i in 1) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  barplot(unlist(ordResult$dev.expl), cex.names=0.8, las=2, ylab="% deviance explained")
}

for (i in 1) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  barplot(unlist(ordResult$RMSE), cex.names=0.8, las=2, ylab="RMSE")
}

#save list of results
phase_space <- data.frame(cbind(ordResult$fits, ordResult$sand_core), ordResult$modspace_fit$Sand,
                          ordResult$modspace_fit$Salinity, ordResult$modspace_fit$Water.depth)
colnames(phase_space) <- c("PC1", "PC2", "traj", "depth", "age", "sand_core", "sand_pred", "salinity_pred",
                           "wdepth_pred")

reg2 <- lm(PC1~sand_core, data=phase_space) 
with(phase_space,plot(PC1, sand_core))
abline(reg2)

reg3 <- lm(PC1~salinity_pred, data=phase_space) 
with(phase_space,plot(PC1, salinity_pred))
abline(reg3)

reg4 <- lm(PC1~wdepth_pred, data=phase_space) 
with(phase_space,plot(PC1, wdepth_pred))
abline(reg4)

plot(sqrt(sand_core),sqrt(ordResult$modspace_fit$Sand), 
     main = "Sand phase space")

abline(lm(sqrt(sand_core)~sqrt(ordResult$modspace_fit$Sand)))

