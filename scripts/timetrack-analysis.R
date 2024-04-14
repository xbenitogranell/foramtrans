##set WD
setwd("/Volumes/xbenitogranell-data/0_project/data/training set")
#setwd("/nfs/xbenitogranell-data/0_project/data/training set")

#clear workspace
rm(list=ls(all=TRUE))
dev.off()
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

##loading libraries for functions used
library(analogue) #to perform timetrack analysis
library(tidyverse) #to manipulate dataframes
library(rioja) #to perform constrained hierarchical clustering
library(mgcv)
library(cluster)
library(ggplot2)
library(ggpalaeo)
library(viridis)
library(ggvegan)
source("/Volumes/xbenitogranell-data/0_project/R codes/useful scripts/GAMs/bennion-frontiers-2015-master/functions.R")
source("/Volumes/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")

#source("/nfs/xbenitogranell-data/0_project/R codes/useful scripts/GAMs/bennion-frontiers-2015-master/functions.R")

## This chunk is to process absolute counts core data

     #read diatom core datasets
      #mergedCores <- read.csv("mergedCores_counts3.csv") #read dataframe with diatom absolute counts including Fondococha
      mergedCores <- read.csv("mergedCores_counts4.csv")[,-1] #read dataframe with diatom absolute counts including Fondococha
      agedepth <- mergedCores[, names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
      diat <- mergedCores[, !names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
      diat[is.na(diat)] <- 0


      diatoms_save <- cbind(agedepth, diat)
      coresList <- split(diatoms_save, diatoms_save$lake)


      # this is function to calculate relative abundance from counts data
      RA <- function(i, cores, ...) {
          core <- cores[[i]]
          core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake"))] # drop year & depths vars
          core <- tran(core, method="percent")
          core[is.na(core)] <- 0
          depth <- coresList[[i]]$depth
          upper_age <- coresList[[i]]$upper_age
          lower_age <- coresList[[i]]$lower_age
          lake <- coresList[[i]]$lake
          cbind.data.frame(depth,upper_age, lower_age, lake, core) #combine extracted columns and remove first row to match with scd

        }

       ## apply RA function function to each core
          coresRA <- lapply(seq_along(coresList), RA, cores=coresList)
          names(coresRA) <- names(coresList)

          #extract dataframes from list
          merged <- plyr::ldply(coresRA, data.frame)
          merged <- merged[,-1] #remove .id variable


          #Select most common species
          agedepth <- merged[, names(merged) %in% c("depth", "upper_age", "lower_age", "lake")]
          diat <- merged[, !names(merged) %in% c("depth", "upper_age", "lower_age", "lake")]

          abund <- apply(diat, 2, max)
          n.occur <- apply(diat>0, 2, sum)
          diat_red <- diat[, n.occur >2 & abund>3] #more than 3% of RA and present in >2 samples

          #check rows with NA
          row.has.na <- apply(diat_red, 1, function(x){any(is.na(x))})
          sum(row.has.na)

          diatoms_save <- cbind(agedepth, diat_red)

          changes <- read.csv("old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
          #new1: ecological groups
          #new2: harmonized taxonomic names

          #transform dataframe to tidy format
          new <- diatoms_save %>%
            gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
            mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
            dplyr::group_by(depth, taxa, lake, upper_age, lower_age) %>%
            summarise(count = sum(count)) %>%
            filter(!count == 0) %>% #this is to remove empty samples (rows)
            filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
            spread(key = taxa, value = count) %>%
            as.data.frame()

          #this is to cutt off to the common era
          cores_merged <- new 
            #mutate(AgeCE = upper_age*(-1)+1950) 
            #filter(AgeCE >= 0)
          
          ## split cores by lakes and reassemble
          coresList <- split(cores_merged, cores_merged$lake)
          
          #save(coresList, file="coresList.RData")

####################################

## This chunk is to process diatom training set

  ##Import data
  environmental_data_lakes <- read.csv("environmental_data_lakes.csv") %>%
            mutate(lake_depth_ratio=Lake_area/Depth_avg) %>%
            mutate(lake_catch_ratio=Lake_area/Wshd_area) %>%
            mutate(catch_vol_ratio=Wshd_area/Vol_total)

        # select.lakes <- paste(c("llaviucu", "yahuarcocha", "pinan", "fondoco"), collapse = '|')
        #   
        #   new <- read.csv("environmental_data_lakes.csv") %>%
        #     filter(str_detect(code, select.lakes))
        #     

  rownames(environmental_data_lakes) <- environmental_data_lakes$code

  training <- read.csv("diatomsTrainingSetDEF.csv", row.names = 1) #with updated diatom taxonomy and selected spp (>3% of RA and present in >2 samples) plus Miriam's Llaviucu slides

  #Regions
  lake_regions <- read.csv("regions.csv", row.names = 1)

  ##Merge training set and regions datasets
  modern_lakes <- merge(training, lake_regions, by="row.names")

  #transform dataframe to tidy format
  df_thin <- modern_lakes %>%
    gather(key = taxa, value = count, -Row.names, -region)#don't gather region

  #import dataframe wiht old and new names to group
  changes_training <- read.csv("old_new_nms_trainingset.csv", stringsAsFactors = FALSE)
  view(changes_training)
  
  #spread
  new <- df_thin %>%
    mutate(taxa = plyr::mapvalues(taxa, from = changes_training$old, to = changes_training$new_1)) %>%
    group_by(region, Row.names, taxa) %>%
    summarise(count = sum(count)) %>%
    spread(key = taxa, value = count)

       
  levels(new$region)

  ##Do some filtering in the training set (i.e., remove datapoints that have spp with way too much abundance)
  #remove.regions <- paste(c("Chile", "Lipez"), collapse = '|')
  
  #this is to reduce training set to northern Andean lakes
  select.regions <- paste(c("Ecuador", "Colombia", "Junin", "Cusco", "eastern"), collapse = '|')
  
  # remove also Titicaca lake as it seems to have a big effect in the ordination scores
  # remove.regions <- paste(c("Chile", "Lipez", "Titicaca"), collapse = '|')
  # remove.lakes <- paste(c("EpNGEO-F_Cubilche1", "Bush-Gd_Miski", "JunnPln_L.Pc,12d1",
  #                         "Cusco_C-PLS-12", "GslngBA_Estrelln", "Brdb-SA_SlcrTrqm", "Bush-Gd_Huamnmrc",
  #                         "GslngBA_KK(nrth)", "Brdb-SA_LgnDsgdr", "EpNGEO-F_Cunrro1"), collapse = '|')
  # 
  # comment spp when analyzing ecological groups data
  new <- new %>%
    #filter(!str_detect(region, remove.regions)) %>%
    #filter(!str_detect(Row.names, remove.lakes)) %>%
    filter(str_detect(region, select.regions)) %>%
    filter(Fragilaria.crotonensis < 40) %>%
    filter(Cyclostephanos.tholiformis < 40) %>%
    filter(Cyclostephanos.andinus < 40) %>%
    as.data.frame()

  training <- new[, -which(names(new) %in% c("Row.names", "region"))]

  #For surface plotting ordination: Merge diatom training set and environmental data of lakes
  row.names(training) <- new[, which(names(new) %in% c("Row.names"))]
  env_surf <- merge(training,environmental_data_lakes, by="row.names")

      # ggplot(env_surf, aes(x=Depth_avg, y=P_B2)) +
      #   geom_smooth(method=lm, se=TRUE)+
      #   geom_point() +
      #   theme_classic()
    
  #For extracting spp from trainingset
  training2 <- env_surf[,2:235]

    #this is for the reduced training set of ecological groups
    #training2 <- env_surf[,2:7]
  
  rowSums(training2)
  
  #For extracting environmental variables from diatom training set
  env_data_lakes <- env_surf[,236:ncol(env_surf)]
  
    #this is for the reduced training set of ecological groups
    #env_data_lakes <- env_surf[,8:ncol(env_surf)]
    
  row.names(env_data_lakes) <- env_data_lakes$code
  
  # For merging environmental dataset with lake regions
  env_data_lakes <- merge(env_data_lakes,lake_regions, by="row.names")
  row.names(env_data_lakes) <- env_data_lakes$code
  env_data_lakes <- env_data_lakes[,-c(1,2)]
  
  
  ## Boxplots
  variables <- c("pH", "Cond", "Water.T", "TP", "Depth_avg", "Ca", "Mg", "K", "Elevation",
                 "MAT", "P.season", "MAP", "T.season", 
                 "Depth_avg", "area_waterbody", "Wshd_area", 
                 "lake_depth_ratio", "lake_catch_ratio", "catch_vol_ratio",
                 "HFP2009","Agriculture", "Crops.and.town", "Grassland.and.shrubs", "region")
  env_data_lakes <- env_data_lakes[,variables]
  
  lakes_wide <- env_data_lakes %>% gather(key=variable, value=value, -region)
  p2 <- ggplot(lakes_wide, aes(x=variable, y=value, fill=region)) + 
    geom_boxplot() +
    facet_wrap(~variable, scale="free")
    
  
  
  
  
  # Comparing Paramo and inter Andean lakes
  ecuador_lakes <- env_data_lakes %>% filter(region %in% c("Ecuador-Andean", "Ecuador-Interandean"))
  ecuador_lakes <- ecuador_lakes[,variables]
  
  ecuador_lakes_wide <- ecuador_lakes %>% gather(key=variable, value=value, -region)
  p2 <- ggplot(ecuador_lakes_wide, aes(x=variable, y=value, fill=region)) + 
    geom_boxplot() +
    facet_wrap(~variable, scale="free")
  
  
   
  
  
      ## Hierarchical Cluster analyses
      diatoms <- training2[, colSums(training2) > 0] #select only present species
      
      #diatoms <- training
      
      # Hellinger transformation
      diatHel <- decostand(diatoms, method="hellinger")
      
      diss <- vegdist(diatHel, method = "bray")
      #clust <- hclust(diss)
      
      pred.flex <- agnes(diss,method='flexible',par.method=c(0.625,0.625,-0.25))
      predflex.hcl <- as.hclust(pred.flex)
      clust <- predflex.hcl
      
      #clust <- chclust(diss, method="coniss")
      #plot(clust,hang=-1, main="Cluster Analysis", cex=0.4)
      
      # Estimate broken stick model
      n <- ncol(diatHel)
      k <- seq(1,n)
      sumFrac <- 1/k
      bstick <- rep(NA,n)
      for(j in 1:n) bstick[j] = 1/n*sum(sumFrac[j:n])
      
      
      # Compare variance explained by each split
      clustVarEx <- rev(diff(clust$height)/max(clust$height))
      plot(k[1:10], clustVarEx[1:10], type="o", col="black", xlab="Number of Splits", ylab="% Var Expl", main="Broken Stick Model")
      points(k[1:10], bstick[1:10], type="o", col="red")
      
      # Split the groups of the ccluster analysis accordingly
      sigClustRef <- cutree(clust, k=4)
      sigClustdf <- data.frame(sigClustRef)
      
      lonLatClust <- cbind(sigClustdf,env_data_lakes) %>%
        select(lat, long, sigClustRef)
      
      

    #This is for prepare data to map diatom community cluster composition
    # FIRST DO CLUSTER ANALYSIS OF DIATOM SPECIES COMPOSITION
      #this is to group species by ecological groups
      new <- df_thin %>% 
        mutate(taxa = plyr::mapvalues(taxa, from = changes_training$old, to = changes_training$new_2)) %>% #ecological grouping
        group_by(region, Row.names, taxa) %>%
        summarise(count = sum(count)) %>%
        filter(!count == "0" ) %>% #this is to remove empty samples (rows)
        ungroup() %>%
        group_by(Row.names, region) %>%
        mutate(relative_abundance_percent = count / sum(count) * 100) %>%
        mutate(plank=sum(count[taxa=="freshwater_planktic" | taxa=="tycoplanktonic"])) %>%
        mutate(benthic=sum(count[taxa=="epiphytics"| taxa== "saline" | taxa=="benthic"])) %>%
        mutate(P_B=plank/benthic) %>%
        mutate(P_B2=(plank-benthic)/(plank+benthic)) %>% #[-1(benthic dominated) to 1(planktic dominated)]
        ungroup() 
      
      #make it wide
      lake_diatom_ratios <- new %>%
        select(region, Row.names, taxa, P_B2, relative_abundance_percent) %>%
        spread(key = taxa, value = relative_abundance_percent) 
      
      lake_diatom_ratios[is.na(lake_diatom_ratios)] <- 0
      
      new <- lake_diatom_ratios
      
      #this is to reduce training set to northern Andean lakes
      select.regions <- paste(c("Ecuador", "Colombia", "Junin", "Cusco", "eastern"), collapse = '|')
      remove.lakes <- paste(c("EpNGEO-F_Cubilche1", "Bush-Gd_Miski", "JunnPln_L.Pc,12d1",
                              "Cusco_C-PLS-12", "GslngBA_Estrelln", "Brdb-SA_SlcrTrqm", "Bush-Gd_Huamnmrc",
                              "GslngBA_KK(nrth)", "Brdb-SA_LgnDsgdr", "EpNGEO-F_Cunrro1"), collapse = '|')

      
      # comment spp when analyzing ecological groups data
      new <- new %>%
        filter(str_detect(region, select.regions)) %>%
        filter(!str_detect(Row.names, remove.lakes)) %>%
        as.data.frame()
      
      training <- new[, -which(names(new) %in% c("Row.names", "region"))]
      
      #For surface plotting ordination: Merge diatom training set and environmental data of lakes
      row.names(training) <- new[, which(names(new) %in% c("Row.names"))]
      env_surf <- merge(training,environmental_data_lakes, by="row.names")
      
      
      row.names(env_surf) <- env_surf[, which(names(env_surf) %in% c("Row.names"))]
      cluster_data <- merge(env_surf, lonLatClust, by="row.names")
      cluster_data <- cluster_data[,-1]
      cluster_data <- cluster_data[,!names(cluster_data) %in% c("lat.y", "long.y", "code")]
    
      top_spp <- cluster_data[,colnames(cluster_data[,2:8])]
      top_spp <- colnames(top_spp)
      top_spp <- match(top_spp, names(cluster_data))
    
      #### Building aggregate composition data for each cluster ####
      ## Averaging over each cluster, and setting each row to sum to one,then taking  column means
      ## code from Pedersen et al 2017
      library(plyr)
      cluster_composition = ddply(cluster_data, .(sigClustRef),
                                  function(x){
                                    x[,2:8] = decostand(x[,2:8],method ="total") 
                                    out_data = data.frame(species = c(names(x)[top_spp]),
                                                          proportion = rep(0,times=length(top_spp)))
                                    for(i in 1:length(top_spp)){
                                      out_data[i,2] = mean(x[,top_spp[i]])
                                    }
                                    return(out_data)
                                  })
    
    cluster_composition <- cluster_composition %>% filter(!species=="P_B2")
  
    cluster_order <- plyr::ddply(cluster_composition,.(sigClustRef), function(x){
      proportion <- sum(x$proportion)
      return(data.frame(proportion=proportion))
    })
    cluster_order = cluster_order$sigClustRef[order(cluster_order$proportion,
                                              decreasing = T)]
    cluster_composition$sigClustRef =factor(cluster_composition$sigClustRef,
                                      levels = cluster_order)
  
  
  # Plot
  cluster_palette <- viridis(4)
    
  #cluster_palette = c("#FFCF4C", "#13ABDA", "#384A81", "#34A576", "#D17634") 
  
  comp_plot <- ggplot(aes(x=sigClustRef, y=proportion),
                     data=cluster_composition)+
    geom_bar(stat="identity", aes(fill=species,order= species))+
    scale_fill_viridis(discrete = TRUE, option = "D") +
    annotate(x=factor(1:4),y = rep(-0.1,times=4),
             colour=cluster_palette,
             geom="point",size=10,shape=18)+
    coord_cartesian(ylim=c(-0.2,1))+
    scale_y_continuous("Average % of community",breaks = c(0,0.25,0.5,0.75,1))+
    scale_x_discrete("Community cluster")+
    theme_bw()+
    theme(text= element_text(size=15),panel.grid.major.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x  = element_blank(),legend.direction="horizontal",
          legend.position = c(0.5,0.9))
  
  #plot modern lake database with clusters of lakes
  library(maps)
  world <- map_data("world")
  southamerica <- ggplot() +
    geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="lightgrey") +
    theme(legend.position = "right")+
    coord_map("albers", parameters = c(-100, -100),  ylim=c(-40,15), xlim=c(-82,-40)) +
    #coord_map("albers", parameters = c(-100, -100),  ylim=c(-10,15), xlim=c(-82,-60)) +
    
    xlab("Longitude") + ylab("Latitude") +
    theme_bw()
  
  legend_title <- "Community cluster"
  
  cluster_map_plt <- southamerica + geom_point(data=lonLatClust, aes(x=long, y=lat, col=factor(sigClustRef)),
                                               shape=18, size=4)+
    scale_color_manual(values = cluster_palette) + 
    scale_fill_discrete(name="Community Cluster") +
    theme(text= element_text(size=15),
          legend.direction="vertical",
          legend.position = "right")
  
  library(cowplot)
  plt <- plot_grid(comp_plot, cluster_map_plt, align="v", axis="tb", nrow = 2,
                   rel_widths = c(1,2))
  
  ggsave("community_clusters_northAndes.png", plt, height = 8, width = 10)
  

####---Do some data exploratory analysis of trainingset

      diatoms <- training2[, colSums(training2) > 0] #select only present species
      
      #diatoms <- training

      # Hellinger transformation
      diatHel <- decostand(diatoms, method="hellinger")
      
      # Run PCA
      diat.pca <- rda(diatHel, scale=T, na.action= "na.exclude", scaling=3)
      # plot(diat.pca, type = 'n', display = 'sites')
      # points(diat.pca, display='sites', col = sigClust, pch = 16)
      # ordihull(diat.pca, groups=sigClust)
      # ordispider(diat.pca, groups = sigClust, label = TRUE)
      # 
      
      # Extract importance of different axes
      pcImp <- summary(diat.pca)$cont$importance

      # Set up the plotting area
      par(mfrow=c(1,1), mar= c(5,4,4,2))

      # select the 15 highest species on each PCA according to the absolute loadings score
      loadPC1 <- order(abs(summary(diat.pca)$species[,1]), decreasing=T)[1:15]
      loadPC2 <- order(abs(summary(diat.pca)$species[,2]), decreasing=T)[1:15]

      # plot the PCA
      plot(summary(diat.pca)$sites[,1], summary(diat.pca)$sites[,2], cex=0.6, xlab="PC1", ylab="PC2", type="n", ylim=c(-1.5,1.4), xlim = c(-1.5,1.5), main= "PCA (Diatoms)" )
      abline(h=0, v=0, lty=2)
      points(summary(diat.pca)$sites[,1], summary(diat.pca)$sites[,2], cex=0.6, xlab="PC1", ylab="PC2", col="dark grey", pch=20 )

      #text(summary(diat.pca)$sites[,1], summary(diat.pca)$sites[,2], labels = rownames(diatoms), cex=0.6)

      specPC1 = summary(diat.pca)$species[,1]
      specPC2 = summary(diat.pca)$species[,2]
      spec = c(loadPC1, loadPC2)
      arrows(x0=rep(0,24), y0=rep(0,24), x1=c(specPC1[spec]), y1=c(specPC2[spec]), col="dark red", length=0.07 )
      text(x=c(specPC1[spec]), y=c(specPC2[spec]), labels=names(specPC1[spec]), cex=0.6)

      # create dataframe
      species_factor_scrs <- data_frame(PCA1=specPC1[spec], PCA2=specPC2[spec])
      species_factor_scrs$species <- names(specPC1[spec])

      
      #Perform NMDS with hellinger species abundance transformation
      diat.nmds <- metaMDS(diatHel, distance="bray", trymax=50, autotransform=F)

      #Plot NMDS
      scrs <- scores(diat.nmds, display = "sites", choices = 1:2)
      plot(diat.nmds, type = "n", xlim=range(scrs[,1]), ylim=range(scrs[,2]))
      
      points(scores(diat.nmds, display="sites"), pch=20)
      
      ordispider(diat.nmds, sigClustRef, conf=0.95, lwd=2,
                        spiders = "centroid", border=1:7, col = 1:7, alpha=63)

      ordihull(diat.nmds, sigClustRef, kind="se", conf=0.95, lwd=2,
              draw = "polygon", col=1:8, border=1:8, alpha=63)

      # ordicluster(diat.nmds, predflex.hcl, prune = 3, col = cutree(predflex.hcl, 8),
      #             draw = "segments")
      # 
      
      #Plot site labels
      plot(diat.nmds, type = "n")
      labels <- as.character(rownames(new))
      text(diat.nmds, labels = labels, pos = 3, cex = 0.4, offset = 0.2)
   

#combine cores and training set (join's analogue package)
df <- analogue::join(coresList[[1]], coresList[[2]], coresList[[3]], coresList[[4]],
                     coresList[[5]], coresList[[6]], coresList[[7]], coresList[[8]], training2, verbose = TRUE)

#check NA in the list
listnans <- lapply(df, function(x) sum(is.na(x)))


#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
  remove <- function(i, cores, ...) {
    core <- cores[[i]]
    core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
    core <- core[, colSums(core) > 0] #select only present species
    # core <- tran(core, "hellinger")
    return(core)
  }

  cores <- lapply(seq_along(df), remove, cores=df)
  #name list elements
  names(cores) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

  #extract training set
  training <- cores$trainingset

  #drop trainingset from the core list
  cores$trainingset <- NULL

  #check NA in the list
  listnans <- lapply(cores, function(x) sum(is.na(x)))


##Perform MULTIPLE time tracks
#tran <- c("sqrt")
tran <- c("hellinger")
#tran <- NULL

doTimeTrack <- function(core, train, scale = TRUE, scaling = 3) {
  timetrack(train, core, method = "rda", scale = scale, transform = tran,
            scaling = scaling, na.rm=TRUE)
}

ttracks <- lapply(cores, doTimeTrack, train = training,
                  scale = FALSE, scaling = 3)


# analogue quality
    #Richard Telford's code
    rlen<-residLen(spp, env, fos, method="cca")
    plot(rlen) # distribution of modern and fossil residual lengths
    plot(depths, rlen$passive, ylab="Squared residual length", xlab="Depth")
    abline(h=quantile(rlen$train, probs = c(0.9,0.95)), col=c("orange", "red"))
    
    
      
    rlen1<-residLen(training_resid, env_data$pH, cores[[lake]], method = "cca")
    rlen2<-residLen(training_resid, env_data$Elevation, cores[[lake]], method = "cca")
    
    depths <- as.numeric(coresList[[lakedepth]]$depth)
    
    rlen1
    rlen2
    plot(rlen) # distribution of modern and fossil residual lengths
    plot(depths, rlen1$passive, ylab="Squared residual length", xlab="Depth")
    plot(depths, rlen2$passive, ylab="Squared residual length", xlab="Depth")
    
    abline(h=quantile(rlen1$train, probs = c(0.9,0.95)), col=c("orange", "red"))
    abline(h=quantile(rlen2$train, probs = c(0.9,0.95)), col=c("orange", "red"))
    

    #merge deleted env data with trainingtset - just once
    training3 <- merge(env_data, training, by="row.names")
    training4 <- training3[,-c(1:18)]
    training_resid <- training4[,colSums(training4) > 0] 
    training_resid <- tran(training_resid, method="hellinger") #give better results transforming
    
    
    par(mfrow=c(4,5))
    par(mar=c(3,3,4,3))
    
    # residual lenghts
    rlenResult <- list()
    
    lake <- "Llaviucu"
    lakedepth <- "llaviucu"
    
    env_data <- env_data_red2
    for(i in 1:length(env_data)){
      rlen<- residLen(training, env_data[,i], cores[[lake]], method = "cca")
      
      rlenResult$passive[[i]] <- rlen$passive
      rlenResult$train[[i]] <- rlen$train
      
      depths <- as.numeric(coresList[[lakedepth]]$depth)
      ages <- as.numeric(coresList[[lakedepth]]$upper_age)
      plot(ages,rlenResult$passive[[i]], 
           ylab="Squared residual length", xlab="Age", main = colnames(env_data[i]))
      
      abline(h=quantile(rlenResult$train[[i]], probs = c(0.9,0.95)), col=c("orange", "red"))
      
    }
    
    names(rlenResult$passive) <- colnames(env_data)
    names(rlenResult$train) <- colnames(env_data)
    
   
    
### CCA
env_data <- env_data_lakes

variables <- c("pH", "Cond", "Water.T", "TP", "Depth_avg", "Ca", "Mg", "K", "Elevation",
               "MAT", "P.season", "MAP", "T.season", 
               "lake_depth_ratio", "lake_catch_ratio", "catch_vol_ratio",
               "HFP2009")

env_data <- env_data[,variables]

#create vector for later use in MAT
lake_depth <- env_data$Depth_avg

#transform variables to meet assumptions of homogenity of variances
env_data <- transform(env_data, Water.T=log10(Water.T+0.25), Elevation=sqrt(Elevation),Cond=log10(Cond+0.25), Ca=log10(Ca+0.25), Mg=log10(Mg+0.25), K=log10(K+0.25), TP=log10(TP+0.25),  
                             Depth_avg=log10(Depth_avg+0.25), MAT=log10(MAT+0.25), P.season=log10(P.season+0.25), MAP=log10(MAP+0.25), T.season=log10(T.season+0.25),
                             lake_depth_ratio=log10(lake_depth_ratio+0.25), lake_catch_ratio=log10(lake_catch_ratio+0.25), catch_vol_ratio=log10(catch_vol_ratio+0.25),
                             HFP2009=log10(HFP2009+0.25))

lake_depth_tran <- env_data$Depth_avg

#plot first cca to subset environmental variables
training <- training[,colSums(training) > 0] 

ccaResult <- list()
for (i in 1:length(env_data)) {
  mod <- cca(training~env_data[,i], na=na.omit, subset = complete.cases(training), scale=TRUE)
  ccaResult$mod[[i]] <- mod
  # plot(ccaResult$mod[[i]], main=colnames(env_data[i]))
  print(anova(ccaResult$mod[[i]]))
}

names(ccaResult$mod) <- colnames(env_data)

#multivariate cca
select <- c("pH", "Cond", "Water.T", "TP", "Ca", "Mg", "K", "Elevation",
               "MAT", "P.season", "MAP", "T.season", 
               "HFP2009")

env_data_red <- env_data[,select]

##boxplots environmental data selected for CCA
par(mfrow = c(4, 4))
par(mar = c(2.5, 3.5, 1, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

for (i in 1:length(env_data_red)) {
  boxplot(env_data_red[,i], cex = 0.6, cex.axis = 0.8,
          las = 1, pch = 19, main=colnames(env_data_red[i])) 
  }

#env_data_red <- env_transformed[,select]
#mod1 <- rda(training~., data=env_data_red, na=na.omit, subset = complete.cases(env_data_red), scale=TRUE)

#Impute missing values with mice package
# CCA env dataset
library(mice)
data <- env_data_red
tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
completedData <- complete(tempData,1)

env_data_red <- completedData
  TP_tran_comp <- env_data_red$TP
  
library(usdm)
vifstep(env_data_red, th=5)

##2nd CCA constrained to environmental data with core samples added passively
#
select <- c("pH", "Cond", "Water.T", "TP", "Ca",
            "MAT", "P.season", "MAP", "T.season", 
            "HFP2009")

env_data_red2 <- env_data_red[,select]

#redo CCA
training <- training[,colSums(training) > 0] 
  
  #duplicate training dataframe before hellinger transformation for subset more abundant species
  training_plt <- training

training <- tran(training, method="hellinger") #give better results transforming
mod_training <- cca(training~., data=env_data_red2, scale=TRUE)
plot(mod_training, scaling = 3)


#Plot eigenvalues and percentages of variation of an ordination object
anova(mod_training, by="axis")
ev <- as.vector(eigenvals(mod_training, model = "constrained")) #extract eigenvalues for then broken stick

evplot <- function(ev) {
  # Broken stick model (MacArthur 1957) Author: Francois Gillet, 25 August 2012
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

br <- evplot(ev)

#
axis.expl <- function(mod, axes = 1:2) {
    if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
}
(labs <- axis.expl(mod_training))
labs[1] <- c(18.1)
labs[2] <- c(16.2)


#make core sample prediction to add passively into CCA
lake <- "Llaviucu"
lakedepth <- "llaviucu"
pred1<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred1 <- as.data.frame(pred1[,1:2])
pred1$lake <- lake
pred1$years <- coresList[[lakedepth]]$upper_age

lake <- "Yahuarcocha"
lakedepth <- "yahuarcocha"
pred2<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred2 <- as.data.frame(pred2[,1:2])
pred2$lake <- lake
pred2$years <- coresList[[lakedepth]]$upper_age

lake <- "Pinan"
lakedepth <- "pinan"
pred3<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred3 <- as.data.frame(pred3[,1:2])
pred3$lake <- lake
pred3$years <- coresList[[lakedepth]]$upper_age

lake <- "Fondococha"
lakedepth <- "fondococha"
pred4<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred4 <- as.data.frame(pred4[,1:2])
pred4$lake <- lake
pred4$years <- coresList[[lakedepth]]$upper_age

pred_cores <- rbind(pred1, pred2, pred3, pred4)
pred_cores <- split(pred_cores, pred_cores$lake)

#png("plt.png", width=10, height=8, units="in", res=300)

# par(mfrow=c(2,2))
# par(mar=c(3,2,2,3))
# 
# for (i in seq_along(pred_cores)) {
#   plot(mod_training, display=c('bp', 'sites'), scaling=3, main=names(pred_cores[i]))
#   points(pred_cores[[i]][,1:2], type="l", col="grey")
#   points(pred_cores[[i]][1, ], pch = 24, cex = 1.6, bg="forestgreen",
#          col = "black")
#   points(pred_cores[[i]][nrow(pred_cores[[i]]),], pch = 22, cex = 1.6,
#          col = "black", bg="forestgreen")
#   legend("topright", pch = c(24, 22), bg=c("forestgreen", "forestgreen"),
#          col=c("black","black"),          cex = 0.8,
#          legend = c("Core top","Core bottom"), bty = "n")
#   print(i)
# }


#Plot for paper
library(viridis)
seq_palette <- viridis(4)

png("CCA_timetrack.png", width=10, height=8, units="in", res=300)

layout(matrix(1:2, ncol = 2))
scrs <- scores(mod_training, display = "species", scaling = 3)

  #abbreviate species names
  spnames_abb <- abbreviate(names(training_plt), minlength=10) 
  colnames(training_plt) <- spnames_abb
  rownames(scrs) <- colnames(training_plt)
  
take <- colnames(training_plt) %in% rownames(scrs)
TAXA <- which(colSums(training_plt[,take] > 0) > 15 & (apply(training_plt[,take]^2, 2, max) > 10)) 
plot(mod_training, display="species", scaling=3, type="n", xlab="", ylab="")
title("Species", adj = 0.3, line = 0.2, cex.main=1, font.main=2, xpd=NA)
title(xlab = paste0(names(labs[1]), " (", sprintf("%.1f", labs[1]), "%)"))
title(ylab = paste0(names(labs[2]), " (", sprintf("%.1f", labs[2]), "%)"))
points(mod_training, display="species", cex = 0.5, scaling = 3, pch=20, select = TAXA)
# ordipointlabel(mod_training, display = "species", scaling = 3,
#                select = TAXA, cex = 0.8, add = TRUE)

spnames_abb <- abbreviate(names(training_plt), minlength=10) 
colnames(training_plt) <- spnames_abb
training_plt <- tran(training_plt, method="hellinger") #give better results transforming
mod_training_plt <- cca(training_plt~., data=env_data_red2, scale=TRUE)
ordipointlabel(mod_training_plt, display = "species", scaling = 3,
               select = TAXA, cex = 0.75, add = TRUE)


plot(mod_training, display=c('bp', 'sites'), xlab="", ylab="")
title("Sites", adj = 0.3, line = 0.2, cex.main=1, font.main=2, xpd=NA)

# title(xlab = paste0(names(labs[1]), " (", sprintf("%.1f", labs[1]), "%)"))
# title(ylab = paste0(names(labs[2]), " (", sprintf("%.1f", labs[2]), "%)"))

points(pred_cores$Pinan[,1:2], type="l", col=seq_palette[1])
points(pred_cores$Yahuarcocha[,1:2], type="l", col=seq_palette[2])
points(pred_cores$Fondococha[,1:2], type="l", col=seq_palette[3])
points(pred_cores$Llaviucu[,1:2], type="l", col=seq_palette[4])


points(pred_cores$Pinan[1, ], pch = 24, cex = 1.6, bg=seq_palette[1],
       col = "grey")
points(pred_cores$Yahuarcocha[1, ], pch = 24, cex = 1.6, bg=seq_palette[2],
       col = "grey")
points(pred_cores$Fondococha[1, ], pch = 24, cex = 1.6, bg=seq_palette[3],
       col = "grey")
points(pred_cores$Llaviucu[1, ], pch = 24, cex = 1.6, bg=seq_palette[4],
       col = "grey")

points(pred_cores$Pinan[nrow(pred_cores$Pinan),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[1])
points(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[2])
points(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[3])
points(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[4])

legend(-1.5, -2.5, pch = c(24, 22), bg=c("forestgreen", "forestgreen"),
       col=c("black","black"),          cex = 0.8,
       legend = c("Core top","Core bottom"), bty = "n")


#par(fig = c(0.02, 0.22, 0.72, 0.97), new = TRUE) #top left
par(fig = c(0.42,0.62, 0.72, 0.97), new = TRUE) #mid left for 1:2 species sires CCA plot

par(mar = c(0, 0, 0, 0))
plot(pred_cores$Pinan[,1:2], type = "n", axes = FALSE, main="")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Pinan[,1:2], type="l", col=seq_palette[1])
points(pred_cores$Pinan[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[1])
points(pred_cores$Pinan[nrow(pred_cores$Pinan),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[1])
text(pred_cores$Pinan[1,], labels=pred_cores$Pinan$years[1], cex = 0.8, offset = 0.5, pos = 3, col="black", 
     bg="white")
text(pred_cores$Pinan[nrow(pred_cores$Pinan),], labels=pred_cores$Pinan$years[nrow(pred_cores$Pinan)], 
     cex = 0.8, offset = 0.5, pos = 3, col="black", bg="white")
title("Piñan", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()


#par(fig = c(0.48,0.68, 0.72, 0.97), new = TRUE) #top right
par(fig = c(0.78,0.98, 0.72, 0.97), new = TRUE) #top right for 1:2 species sires CCA plot

par(mar = c(0, 0, 0, 0))
plot(pred_cores$Yahuarcocha[,1:2], type = "n", axes = FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Yahuarcocha[,1:2], type="l", col=seq_palette[2])
points(pred_cores$Yahuarcocha[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[2])
points(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[2])
text(pred_cores$Yahuarcocha[1,], labels=pred_cores$Yahuarcocha$years[1], cex=0.8, offset=0.5, pos=4, col="black")
text(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], labels=pred_cores$Yahuarcocha$years[nrow(pred_cores$Yahuarcocha)],
     cex = 0.8, offset = 0.5, pos = 4, col="black")
title("Yahuarcocha", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

#par(fig = c(0.02, 0.22, 0.02, 0.27), new = TRUE) #bottom left
par(fig = c(0.42,0.62, 0.02, 0.27), new = TRUE) #mid bottom for 1:2 species sires CCA plot
par(mar = c(0, 0, 0, 0))
plot(pred_cores$Fondococha[,1:2], type = "n", axes=FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Fondococha[,1:2], type="l", col=seq_palette[3])
points(pred_cores$Fondococha[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[3])
points(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[3])
text(pred_cores$Fondococha[1,], labels=pred_cores$Fondococha$years[1], cex=0.8, offset=0.5, pos=3, col="black")
text(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], labels=pred_cores$Fondococha$years[nrow(pred_cores$Fondococha)],
     cex=0.8, offset = 0.5, pos = 4, col="black")
title("Fondococha", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

#par(fig = c(0.48,0.68, 0.02, 0.27), new = TRUE) #bottom right
par(fig = c(0.78,0.98, 0.02, 0.27), new = TRUE) #bottom right for 1:2 species sires CCA plot
par(mar = c(0, 0, 0, 0))
plot(pred_cores$Llaviucu[,1:2], type = "n", axes = FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Llaviucu[,1:2], type="l", col=seq_palette[4])
points(pred_cores$Llaviucu[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[4])
points(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[4])
text(pred_cores$Llaviucu[1,], labels=pred_cores$Llaviucu$years[1], cex=0.8, offset=0.5, pos=2, col="black")
text(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], labels=pred_cores$Llaviucu$years[nrow(pred_cores$Llaviucu)], 
     cex=0.8, offset=0.5, pos=3, col="black")
title("Llaviucu", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

dev.off()


#Modern analogue technique
MATResult <- list()

# for(i in 1:length(env_data_red2)){
#   mod <- rioja::MAT(training, env_data_red2[,i])
#   MATResult$mod[[i]] <- mod
#   pred <- predict(MATResult$mod[[i]], cores[[lake]])
# 
#   MATResult$pred[[i]] <- pred #closest analogue
#   
#   depths <- as.numeric(coresList[[lakedepth]]$depth)
#   ages <- as.numeric(coresList[[lakedepth]]$upper_age)
#   
#   plot(ages,MATResult$pred[[i]]$dist.n[,1], 
#        ylab="Squared chord distance", xlab="Age", main = colnames(env_data[i]))
#   
#   goodpoorbad<-quantile(paldist(training), prob=c(0.05, 0.1))
#   abline(h=goodpoorbad, col=c("orange", "red"))
# }



## 
env <- env_data$Depth_avg
lake <- "Fondococha"
lakedepth <- "fondococha"

mod <- rioja::MAT(training/100, env, dist.method="sq.chord")
pred <- predict(mod,cores[[lake]]/100)
ages <- as.numeric(coresList[[lakedepth]]$upper_age)
plot(ages, pred$dist.n[,1], ylab="Squared chord distance", xlab="Age cal yr BP")
goodpoorbad <- quantile(pred$dist.n[,1], probs = c(0.75, 0.95))
abline(h=goodpoorbad, col=c("orange", "red"))

dist_to_analogues <- as.data.frame(ages) %>% 
  mutate(
    dist_to_analogues = pred$dist.n[, 1],
    quality  = cut(dist_to_analogues, breaks = c(0, goodpoorbad, Inf), labels = c("good", "poor", "bad"))
  )
attr(dist_to_analogues, which = "goodpoorbad") <- goodpoorbad

dist_to_analogues_plot <- ggpalaeo:::plot_diagnostics(x = dist_to_analogues, x_axis = "ages", y_axis = "dist_to_analogues", 
                                                      goodpoorbad = attr(dist_to_analogues, "goodpoorbad"), fill = c("salmon", "lightyellow", "skyblue"), categories = c("Good", "Fair", "None")) + 
  labs(x = "Cal yr BP", y = "Squared chord distance", fill = "Analogue quality") +
  ggtitle("Fondococha")

ggsave("dist_to_analogues_plot_fondococha.png", dist_to_analogues_plot, height = 8, width = 10)

# dist_to_analogues$lake <- "Pinan"
# dist_to_analogues_pinan <- dist_to_analogues
# 
# dist_to_analogues$lake <- "Yahuarcocha"
# dist_to_analogues_yahuarcocha <- dist_to_analogues
# 
# dist_to_analogues$lake <- "Fondococha"
# dist_to_analogues_fondococha <- dist_to_analogues
# 
# dist_to_analogues$lake <- "Llaviucu"
# dist_to_analogues_llaviucu <- dist_to_analogues

dist_to_analogues_all <- bind_rows(dist_to_analogues_pinan, dist_to_analogues_yahuarcocha,
                                   dist_to_analogues_fondococha, dist_to_analogues_llaviucu)

#rename Piñan
dist_to_analogues_all <- dist_to_analogues_all %>% mutate(lake=str_replace(lake,"Pinan", "Piñan"))

#arrange dataframe by latitude of lakes
dist_to_analogues_all$lake <- factor(dist_to_analogues_all$lake, levels = c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu"))

write.csv(dist_to_analogues_all, "dist_to_analogues_all_lakes.csv")

dist_to_analogues_all_plt <- ggplot(data=dist_to_analogues_all, aes(x=ages, y=dist_to_analogues, col=quality)) +
  geom_point() +
  scale_color_manual(values=c("skyblue", "orange", "red"))+
  facet_wrap(~lake, scales = "free")+
  xlab("Cal yr BP") +
  ylab("Squared chord distance") +
  labs(col = "Quality")+
  theme_bw()

ggsave("dist_to_analogues_plot_all.png", dist_to_analogues_all_plt, height = 8, width = 10)

##


#plot using ggvegan package
autoplot(mod1)

fdat <- fortify(mod1) 

ford <- fortify(mod1, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

mul2 <- ggvegan:::arrowMul(arrows[, take], pred)


#plot
ggplot() +
  geom_point(data = subset(ford, Score == 'sites'),
             mapping = aes(x = CCA1, y = CCA2)) + 
  geom_segment(data = arrows,
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1)) +
  #geom_line(data=pred*mul2, aes(x=CCA1, y=CCA2), colour="lightblue")+
  coord_fixed()


## plot timetracks
par(mfrow=c(3,3))
par(mar=c(2,2,2,2))

for(i in seq_along(ttracks)) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  plot(ttracks[[i]], type="n")
  points(ttracks[[i]], which = "ordination", pch = 20, col="grey", cex = 0.7)
  #points(ttracks[[i]], which = "ordination", col = sigClustRef, pch = 20, cex = 0.7)
  # with(training, ordihull(ttracks[[i]], which = "ordination", groups = sigClustRef,
  #          draw="polygon", col=1:4, border=1:4, alpha = 63, lwd=1))
  # #ordihull(ttracks[[i]]$ordination, display = "sites", groups=sigClust, draw = "polygon", col=1:4)
  title(names(ttracks[i]))
  par(op)

}

## Plot most common/abundant taxa
layout(matrix(1:2, ncol = 2))

op <- par(mar = c(5,4,1,1) + 0.1)
set.seed(24)
scrs <- scores(ttracks[[1]]$ordination, display = "species", scaling = 3)
  scrs <- scores(mod_training, display = "species", scaling = 3)

take <- colnames(training) %in% rownames(scrs)
TAXA <- which(colSums(training[,take] > 0) > 15 & (apply(training[,take]^2, 2, max) > 10) &
                ((scrs[,1] <= -0.35 | scrs[,1] >= 0.5) | (scrs[,2] <= -0.35 | scrs[,2] >= 0.25)))

plot(ttracks[[1]]$ordination, display = "species", scaling = 3, type = "n")
  plot(mod_training, display="species", scaling=3, type="n")
  
# plot(ttracks[[1]]$ordination, display = "sites", scaling = 3, type = "n")
# scrs <- scores(ttracks[[1]]$ordination, display = "sites", scaling = 3)
# points(ttracks[[1]]$ordination, display = "sites", cex = 0.5, scaling = 3, col=sigClust, pch=20)
# ordispider(ttracks[[1]]$ordination, sigClust, spiders = "centroid")
#  

rect(-0.35, -0.35, 0.5, 0.25, border = "grey", col = NA)
points(ttracks[[1]]$ordination, display = "species", cex = 0.5, scaling = 3)
  #plot for CCA
  points(mod_training, display="species", cex = 0.5, scaling = 3, pch=20, select = TAXA)
  ordipointlabel(mod_training, display = "species", scaling = 3,
                 select = TAXA, cex = 0.5, add = TRUE)
  
ordipointlabel(ttracks[[1]]$ordination, display = "species", scaling = 3,
               select = TAXA, cex = 0.5, add = TRUE)
  

#### Plot rest of the taxa
TAXA <- which(colSums(training[,take] > 0) > 5 & (apply(training[,take]^2, 2, max) > 5))
plot(ttracks[[1]]$ordination, display = "species", scaling = 3, type = "n",
     xlim = c(-0.35, 0.5), ylim = c(-0.35, 0.25))
rect(-0.35, -0.35, 0.5, 0.25, border = "grey", col = NA)
points(ttracks[[1]]$ordination, display = "species", cex = 0.5, scaling = 3)
ordipointlabel(ttracks[[1]]$ordination, display = "species", scaling = 3,
               xlim = c(-1,1.75), ylim = c(-0.75,0.75),
               select = TAXA, cex = 0.5, add = TRUE)

# project biplot vectors
scrs <- scores(ttracks[[1]]$ordination, display = "sites", scaling = 3)
plot(ttracks[[1]]$ordination, display = "sites", scaling = 3, type = "n")
dat <- as.data.frame(with(env_transformed, cbind(scrs, TP, pH, 
      Cond, K, Ca, Depth_avg, lake_depth_ratio, T.season, Elevation, MAT, Water.T, P.season, HFP2009)))
diat.ev <- envfit(scrs ~., data = dat[,3:ncol(dat)], na.rm=TRUE)
plot(diat.ev)


dev.off()


# taxa present in more than 20 lakes with >10% abundance relative
TAXA <- which(colSums(training[,take] > 0) > 20 & (apply(training[,take]^2, 2, max) > 10))
scrs <- scores(ttracks[[1]]$ordination, display = "species", scaling = 3)

# create dataframe
species_factor_scrs <- data_frame(PCA1=scrs[TAXA], PCA2=scrs[TAXA])
scrsdf <- as.data.frame(scrs)
scrsdf$species <- rownames(scrsdf)
rownames(scrsdf) <- NULL

species_factor_scrs <- scrsdf[TAXA,]
      
plot(species_factor_scrs[,1], species_factor_scrs[,2])
text(species_factor_scrs[,1], species_factor_scrs[,2], species_factor_scrs[,3], cex=0.5, offset = 1)

# create voronoi line sigments
voronoi <- deldir(species_factor_scrs$PC1, species_factor_scrs$PC2)

ggplot(data=species_factor_scrs, aes(x=PC1,y=PC2)) +
  #Plot the voronoi lines
  geom_segment(
    aes(x = x1, y = y1, xend = x2, yend = y2),
    size = 2,
    data = voronoi$dirsgs,
    linetype = 1,
    color= "#FFB958") + 
  #Plot the points
  geom_point(
    fill=rgb(70,130,180,255,maxColorValue=255),
    pch=21,
    size = 4,
    color="#333333")


## This is to make a loop for ordisurf variables
variables <- c("pH", "Cond", "Water.T", "TP", "Depth_avg", "Ca", "Mg", "K", "Elevation",
               "MAT", "P.season", "MAP", "T.season", 
               "lake_depth_ratio", "lake_catch_ratio", "catch_vol_ratio",
               "HFP2009")

matrix <- env_surf[,variables]

    #transform variables to meet assumptions of homogenity of variances
    env_transformed <- transform(matrix, Water.T=log10(Water.T+0.25), Elevation=sqrt(Elevation),Cond=log10(Cond+0.25), Ca=log10(Ca+0.25), Mg=log10(Mg+0.25), K=log10(K+0.25), TP=log10(TP+0.25),  
                                 Depth_avg=log10(Depth_avg+0.25), MAT=log10(MAT+0.25), P.season=log10(P.season+0.25), MAP=log10(MAP+0.25), T.season=log10(T.season+0.25),
                                 lake_depth_ratio=log10(lake_depth_ratio+0.25), lake_catch_ratio=log10(lake_catch_ratio+0.25), catch_vol_ratio=log10(catch_vol_ratio+0.25),
                                 HFP2009=log10(HFP2009+0.25))
    
    
    #panel correlation plots to assess data distribution
    panel.hist <- function(x, ...) {     
      usr <- par("usr"); on.exit(par(usr))     
      par(usr = c(usr[1:2], 0, 1.5) )     
      h <- hist(x, plot = FALSE)     
      breaks <- h$breaks; nB <- length(breaks)     
      y <- h$counts; y <- y/max(y)     
      rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
    }
    
    panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
      usr <- par("usr"); on.exit(par(usr))     
      par(usr = c(0, 1, 0, 1))     
      r <- abs(cor(x, y, use = "complete"))   
      txt <- format(c(r, 0.123456789), digits = digits)[1]     
      txt <- paste0(prefix, txt)     
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
      text(0.5, 0.5, txt, cex = cex.cor * r) }
    
    
    # Plot
    pairs(env_transformed, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
    

## This is to resolve the problem when the target env_surf variable has no NAs;
# insert NAs randomly
HFP2009 <- env_surf$HFP2009
ind <- which(HFP2009 %in% sample(HFP2009, 5))
HFP2009[ind]<-NA
env_surf$HFP2009 <- HFP2009

surf <- ordisurf(ttracks[[1]]$ord, env_surf$HFP2009,
                 method = "REML", scaling = ttracks[[1]]$scaling,
                 select = FALSE, plot = TRUE, knots = 20)

summary(surf)



## fit univariate surface - only once
surf <- ordisurf(ttracks[[1]]$ord, env_transformed[,"TP"],
                 method = "REML", scaling = ttracks[[1]]$scaling,
                 select = TRUE, plot = FALSE, knots = 20)


#
RMSE <- function(error) { sqrt(mean(error^2)) }
RMSE(surf$residuals)

# plot timetracks using Gavin's function
par(mfrow=c(2,2))
par(mar=c(3,3,2,3))

#surf <- NULL
tt.lwd <- 1

lakes <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")
ttracks <- ttracks[lakes]


for(i in seq_along(ttracks[lakes])) {
  #png("ecuador_lakes_timetrack.png", width =12, height=6, units="in", res=300)
  ttplot(ttracks[[i]], surf = surf, lwd = tt.lwd)
  title(names(ttracks[i]))
  legend("topleft", col="#3465a4", lty=1, cex=1,
         legend= c("Total Phosphorous"), bty='n')
  
}


##constrained timetrack
tran <- "sqrt"
doTimeTrackCons <- function(core, train, env, formula, scale = TRUE, scaling = 3) {
  timetrack(train, core, method = "cca", scale = scale, transform = tran,
            scaling = scaling, na.rm=TRUE)
}


ttracksCons <- lapply(cores, doTimeTrackCons, env=env_data, formula=~., train = training,
                      scale = TRUE, scaling = 3)


##
tt.lwd <- 1
#surf <- NULL

lake <- c("Llaviucu")
lake_ages <- c("llaviucu")
k=3 #number of groups in ccluster

for (i in seq_along(ttracks[lake])) {
    tpch = 21
    bpch = 22
    pcol = "#a40000"
    pbg = "#ef2929"
    pcex = 1.6
    
    apch = 24
    acol = "orange"
    abg = "orange"
    
    multittplot(ttracks[[lake]], lwd = tt.lwd, surf=surf)
    title(lake)
    legend("topleft", col="#3465a4", lty=1, cex=1,
           legend= c("Total Phosphorous"), bty='n')
    
    coreCluster_lake <- coreCluster[[lake]]
    sigClust <- cutree(coreCluster_lake, k=k)
    agesClust <- as.numeric(coresList[[lake_ages]]$upper_age)
    fits <- as.data.frame(cbind(fitted(ttracks[[lake]], type = "passive"), 
                                traj= sigClust, brks = agesClust))
    
    # surf <- ordisurf(ttracks[[1]]$ord, env_data[,i], method = "REML", scaling = ttracks[[1]]$scaling,
    #                  select = TRUE, add = TRUE, main = colnames(env_data[i]), knots = 20, 
    #                  npoints = nrow(fits), col = "#3c73ad") #npoints is the number of locations at which to evaluate the fitted surface
    # 
    
    segments(fits$PC1[-length(fits$PC1)],
             fits$PC2[-length(fits$PC2)],
             fits$PC1[-1L],fits$PC2[-1L], col=fits$traj, type="b")
    
    locate <- cumsum(rle(fits$traj)$lengths)+1
    points(fits[locate, ], pch=24, col=acol, bg=abg)
    text(fits[locate, ], labels=fits$brks[locate], cex = 1, offset = 0.5, pos = 4)
}



## Perform MULTIPLE constrained cluster analyses
doCluster <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- decostand(core, method="hellinger") #Hellinger transform relative abundance data
  diss <- vegan::vegdist(core, method="bray")
  clust <- chclust(diss, method="coniss")
  #bstick(clust)
  clust
}

coreCluster <- lapply(seq_along(cores), doCluster, cores = cores)
names(coreCluster) <- names(cores)

## plot Clusters
par(mfrow=c(3,3))
par(mar=c(4,2,2,2))

for(i in seq_along(coreCluster)) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  plot(coreCluster[[i]], hang=-1)
  title(names(coreCluster[i]))
  par(op)
}

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
    # a splits = a+1 zones (k), so split the groups of the ccluster analysis accordingly
    sigClustLlav <- cutree(coreCluster$Llaviucu, k=4)
    locate <- cumsum(rle(sigClustLlav)$lengths)+1
    zones_llaviucu <- coresList$llaviucu$upper_age[locate]

    sigClustPin <- cutree(cCluster_Pinan, k=4)
    locate <- cumsum(rle(sigClustPin)$lengths)+1
    zones_pinan <- coresList$pinan$upper_age[locate]

    #sigClustTriumfo <- cutree(cCluster_Triumfo, k=3)
    #sigClustUmy <- cutree(cCluster_Umayo, k=3)
    sigClustYah <- cutree(cCluster_Yahuarcocha, k=4)
    locate <- cumsum(rle(sigClustYah)$lengths)+1
    zones_yahuarcocha <- coresList$yahuarcocha$upper_age[locate]
    
    sigClustFondo <- cutree(cCluster_Fondococha, k=4)
    locate <- cumsum(rle(sigClustFondo)$lengths)+1
    zones_fondococha <- coresList$fondococha$upper_age[locate]
    
    sigClustTiti <- cutree(cCluster_Titicaca, k=3)


# separate cores by regions
    north <- coresList[c(1,3,4,8)]
    central <- coresList[c(2,5,7)]

    north_diat <- cores[c(1,3,4,8)]
    central_diat <- cores[c(2,5,7)]

#this is to check temporal resolution
    par(mfrow=c(3,3))

    dif <- c(NA, diff(north[["llaviucu"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(north[["llaviucu"]][["upper_age"]])
    title("Llaviucu")

    dif <- c(NA, diff(north[["yahuarcocha"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(north[["yahuarcocha"]][["upper_age"]])
    title("Yahuarcocha")

    dif <- c(NA, diff(north[["pinan"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(north[["pinan"]][["upper_age"]])
    title("Pinan")

    dif <- c(NA, diff(north[["fondococha"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(north[["fondococha"]][["upper_age"]])
    title("Fondococha")

    dif <- c(NA, diff(central[["titicaca"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(central[["titicaca"]][["upper_age"]])
    min(central[["titicaca"]][["upper_age"]])
    title("Titicaca")

    dif <- c(NA, diff(central[["umayo"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(central[["umayo"]][["upper_age"]])
    min(central[["umayo"]][["upper_age"]])
    title("Umayo")

    dif <- c(NA, diff(central[["lagunillas"]][["upper_age"]])) #NA
    median(dif[-1])
    plot.ts(dif)
    max(central[["lagunillas"]][["upper_age"]])
    min(central[["lagunillas"]][["upper_age"]])
    title("Lagunillas")


#Binning samples by categories of ages
doBin <- function(i, cores, ...) {
  # Function for binning samples in each core (Alistair Seddon)
  source("/Volumes/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  core <- cores[[i]]
  core <- core[, colSums(core) > 0] #select only present species

  #ages <- north[[i]]$upper_age #here could be upper_age or AgeCE ###NORTH ANDES
  #ages <- round(central[[i]]$upper_age, digits = 0) #here could be upper_age or AgeCE ###CENTRAL ANDES
  ages <- round((north[[i]]$upper_age*(-1)+1950), digits=0) #this is to get ages_AD
  
  #make the age categories
  diatomBin1 <- binFunc(as.data.frame(core), as.numeric(ages), 100, 0, 1700) ### NORTH lakes
  diatomBin2 <- binFunc(as.data.frame(core), as.numeric(ages), 30, 1700, 2000) ##

  #merge the two binning dataframes for North lakes
  rbind.data.frame(diatomBin1, diatomBin2)

  #this is to remove rows with NaNs after binning
  # row.all.na <- apply(diatomBin, 1, function(x) all(is.na(x)))
  # sum(row.all.na)
  # diatomBin <- diatomBin[ !row.all.na, ]
  # diatomBin
 }

coreBinned_north <- lapply(seq_along(north_diat), doBin, cores = north_diat) ###NORTH ANDES
coreBinned_central <- lapply(seq_along(central_diat), doBin, cores = central_diat) ###CENTRAL ANDES

#check NA in the list
listnans <- lapply(coreBinned_north, function(x) sum(is.na(x)))

names(coreBinned_north) <- names(north_diat)
names(coreBinned_central) <- names(central_diat)


#this is to extract binned ages; only needs to do it once because binned ages are the same for all the cores
ages <- as.data.frame(as.numeric(row.names(coreBinned_north$Llaviucu)))
ages <- plyr::rename(ages,c("as.numeric(row.names(coreBinned_north$Llaviucu))"="age"))

ages <- as.data.frame(as.numeric(row.names(coreBinned_central$Lagunillas)))
ages <- plyr::rename(ages,c("as.numeric(row.names(coreBinned_central$Lagunillas))"="age"))



#do interpolation between adjacent samples
library(zoo)
doInterpol <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- as.data.frame(na.approx(core, na.rm = TRUE)) #do interpolation between adjacent samples

  #merge with ages vector
  core <- merge(core, ages, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names"); after interpolation, rownames are re-started
  core <- core[,-1] #remove first column
  rownames(core) <- core$age #put ages as rownames
  core <- core[order(core$age),] #order by increasing age
  core <- core[,!names(core) %in% c("age")] #drop age variable

  #this is to remove rows with NaNs after merging
  # row.has.na <- apply(core, 1, function(x){any(is.na(x))})
  # sum(row.has.na)
  # core <- core[!row.has.na,]
  # core

  #this is to insert zero with NaNs after merging which results in the same number of rows (years)
  #core[is.na(core)] <- 0
  #core
}

coreBinned_north_interPol <- lapply(seq_along(coreBinned_north), doInterpol, cores = coreBinned_north) ###NORTH ANDES
names(coreBinned_north_interPol) <- names(north_diat)

coreBinned_central_interPol <- lapply(seq_along(coreBinned_central), doInterpol, cores = coreBinned_central) ###CENTRAL ANDES
names(coreBinned_central_interPol) <- names(central_diat)


      #extract dataframes from list
      ages <- plyr::ldply(as.numeric(rownames(coreBinned_north_interPol[[1]])), data.frame)
      ages <- plyr::rename(ages,c("X..i.."="age"))

      diat <- plyr::ldply(coreBinned_north_interPol, data.frame)
      diat <- plyr::rename(diat,c(".id"="lake"))
      diat <- cbind(ages,diat)

      write.csv(diat, "diat_binned_data.csv")
      
      #now the dataframe is ready to filter by ages
      diat_filt <- diat %>% filter(age<1100)
      diatList <- split(diat_filt, diat_filt$lake)


#check NA in the list
listnans <- lapply(coreBinned_north_interPol, function(x) sum(is.na(x)))



  #something odd with Lagunillas
  lagunillas <- central_diat$Lagunillas
  lagunillasBin <- binFunc(as.data.frame(lagunillas), central$lagunillas$upper_age, 30, 1500, 7000)

  rownames(fossildata) <- rownames(lagunillasBin[-1,])


#rename list to do timetrack and betadiversity calculations aftwerwards
coreBinned_north <- coreBinned_north_interPol
coreBinned_central <- coreBinned_central_interPol

coreBinned <- coreBinned_north
coreBinned <- coreBinned_central

## do timetrack with core binned data

ttracks_binned <- lapply(coreBinned, doTimeTrack, train = training,
                      scale = FALSE, scaling = 3)

    # plot timetrack binned data
    #surf <- NULL
    tt.lwd <- 1

    for(i in seq_along(ttracks_binned)) {
      ttplot(ttracks_binned[[i]], surf = NULL, lwd = tt.lwd)
      title(names(ttracks_binned[i]))
    }


#do cluster analysis with binned data
coreClusterBinNorth <- lapply(seq_along(coreBinned_north), doCluster, cores = coreBinned_north)
names(coreClusterBinNorth) <- names(coreBinned_north)

coreClusterBinCentral <- lapply(seq_along(coreBinned_central), doCluster, cores = coreBinned_central)
names(coreClusterBinCentral) <- names(coreBinned_central)

#names(coreClusterBin) <- names(cores)

coreClusterBin <- coreClusterBinNorth
coreClusterBin <- coreClusterBinCentral

## plot Clusters
par(mfrow=c(3,3))

for(i in seq_along(coreClusterBin)) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  plot(coreClusterBin[[i]], hang=-1)
  title(names(coreClusterBin[i]))
  par(op)
}


    # Estimate broken stick model for significant zones (visual inspection) on core binned data
    ## plot Clusters
    par(mfrow=c(3,3))
    #par(mar=c(1,1,1,1))

    #png(filename="plot.png", width=11.69, height=8.27, units=in, res=300)


    for (i in seq_along(coreClusterBin)) {
      core <- coreClusterBin[[i]]
      n <- ncol(cores[[i]])
      k <- seq(1,n)
      sumFrac <- 1/k
      bstick <- rep(NA,n)
      for(j in 1:n) bstick[j] = 1/n*sum(sumFrac[j:n])

      # Compare variance explained by each split
      clustVarEx <- rev(diff(core$height)/max(core$height))
      plot(k[1:10], clustVarEx[1:10], type="o", col="black", xlab="Number of Splits", ylab="% Var Expl")
      title(names(coreClusterBin[i]))
      points(k[1:10], bstick[1:10], type="o", col="red")

    }


    dev.off()


    #extract core binned clusters
    nams <- names(coreClusterBin)
    for (i in seq_along(coreClusterBin)) {
      assign(paste0("cClusterBin_", nams[i]), coreClusterBin[[i]])
    }


    # a splits = a+1 zones (k), so split the groups of the ccluster analysis accordingly
    sigClustLlav <- cutree(coreClusterBin$Llaviucu, k=3)
    sigClustPin <- cutree(coreClusterBin$Pinan, k=3)
    sigClustTriumfo <- cutree(coreClusterBin$Triumfo, k=3)
    sigClustUmy <- cutree(coreClusterBin$Umayo, k=4)
    sigClustYah <- cutree(coreClusterBin$Yahuarcocha, k=4)
    sigClustTiti <- cutree(coreClusterBin$Titicaca, k=3)
    sigClustLagu <- cutree(coreClusterBin$Lagunillas, k=3)
    sigClustFondo <- cutree(coreClusterBin$Fondococha, k=3)
    

    #sigClust <- cutree(coreCluster$Lagunillas, k=3)


    # Find the age intervals of the zones identified by the cluster analysis
    age <- row.names(coreBinned$Llaviucu)
    age <- row.names(coreBinned$Yahuarcocha)
    age <- row.names(coreBinned$Pinan)
    age <- row.names(coreBinned$Triumfo)
    age <- row.names(coreBinned$Umayo)
    age <- row.names(coreBinned$Titicaca)
    age <- row.names(coreBinned$Fondococha)

    # name vector
    names(sigClustLlav) <- age
    names(sigClustYah) <- age
    names(sigClustPin) <- age
    names(sigClustTriumfo) <- age
    names(sigClustUmy) <- age
    names(sigClustTiti) <- age
    names(sigClustFondo) <- age
    
    #
    agesClust <- as.numeric(names(sigClustLlav))
    agesClust <- as.numeric(names(sigClustYah))
    agesClust <- as.numeric(names(sigClustPin))
    agesClust <- as.numeric(names(sigClustTriumfo))
    agesClust <- as.numeric(names(sigClustUmy))
    agesClust <- as.numeric(names(sigClustTiti))
    agesClust <- as.numeric(names(sigClustFondo))


## build timetrack plots

    dev.off()
    png("timetrack.png", width = 11, height = 8, res = 300, units = "in")
    par(mfrow=c(2,2))

    par(mar=c(4,4,1,1.5))

    ttrack <- ttracks$Llaviucu
    ttrack <- ttracks_binned$Llaviucu

    ttrack <- ttracks$Yahuarcocha
    ttrack <- ttracks_binned$Yahuarcocha

    ttrack <- ttracks$Pinan
    ttrack <- ttracks_binned$Pinan

    ttrack <- ttracks$Fondococha
    ttrack <- ttracks_binned$Fondococha

    ttrack <- ttracks$Titicaca
    ttrack <- ttracks_binned$Titicaca

    ttrack <- ttracks$Umayo
    ttrack <- ttracks_binned$Umayo

    ttrack <- ttracks$Triumfo
    ttrack <- ttracks_binned$Triumfo

    ttrack <- ttracks$Triumfo
    ttrack <- ttracks_binned$Lagunillas


    #this is for binned and non-binned core data
    sigClust <- sigClustLlav
    sigClust <- sigClustYah
    sigClust <- sigClustPin
    sigClust <- sigClustFondo
    sigClust <- sigClustTiti
    sigClust <- sigClustUmy
    sigClust <- sigClustTriumfo
    sigClust <- sigClustLagu


    # Find the age intervals of the zones identified by the cluster analysis
    # this is for non-binned core data
    agesClust <- as.numeric(coresList$llaviucu$upper_age)
    agesClust <- as.numeric(coresList$pinan$upper_age)
    agesClust <- as.numeric(coresList$yahuarcocha$upper_age)
    agesClust <- as.numeric(coresList$fondococha$upper_age)

    agesClust <- as.numeric(coresList$triumfo$upper_age)
    agesClust <- as.numeric(coresList$umayo$upper_age)
    agesClust <- as.numeric(coresList$titicaca$upper_age)
    agesClust <- as.numeric(coresList$lagunillas$upper_age)

    
    fits <- as.data.frame(cbind(fitted(ttrack, type = "passive"), traj= sigClust, brks = agesClust)) #passive

      #this is to cutt to the first rows (after binning something odd appears in that fits dataframe duplicates)
      #fits <- fits[1:length(age),]

    scrs <- scores(ttrack$ordination, choices = 1:2, scaling = ttrack$scaling, display = "sites")

    xlim <- range(scrs[, 1], fits[, 1])
    ylim <- range(scrs[, 2], fits[, 2])

    plot(ttrack, type = "n", ptype = "n", xlim = xlim, ylim = ylim, main="", xlab = "PCA1", ylab = "PCA2")
    points(scores(ttrack, type = "passive"),col = "grey", pch = 19)

        # project biplot vectors
        dat <- as.data.frame(with(env_transformed, cbind(scrs, TP, pH, Cond, K, Ca, Depth_avg, lake_depth_ratio)))
        diat.ev <- envfit(scrs ~ TP + Ca + pH + Cond + Depth_avg + lake_depth_ratio, data = dat, na.rm=TRUE)
        plot(diat.ev)
        
    #this is to colour reference spatial association of lakes
    #points(scores(ttrack, type = "passive"),col = as.integer(sigClustdf$sigClust), pch = 19)


    surf1 <- ordisurf(ttrack$ord, env_surf$TP,
                   method = "REML", scaling = ttrack$scaling, knots = 20,
                   col = "#3465a4", plot = TRUE, lwd.cl = 1, cex = 3, add = TRUE)
    summary(surf1)
    TP.fitted <- surf1$fitted.values


    #plot training set by lake regions
    #points(scores(ttrack, type = "passive"), col=as.integer(modern_lakes$region), pch = 19)


    segments(fits$PC1[-length(fits$PC1)],
             fits$PC2[-length(fits$PC2)],
             fits$PC1[-1L],fits$PC2[-1L], col=fits$traj, type="b")

    tpch = 21
    bpch = 22
    pcol = "#a40000"
    pbg = "#ef2929"
    pcex = 1.6

    apch = 24
    acol = "orange"
    abg = "orange"

    points(fits[1, ], pch = tpch, cex = pcex,
           col = pcol, bg = pbg)
    points(fits[nrow(fits), ], pch = bpch, cex = pcex,
           col = pcol, bg = pbg)

    #identify where trajectories change in the fits dataframe
    locate <- cumsum(rle(fits$traj)$lengths)+1
    points(fits[locate, ], pch=24, col=acol, bg=abg)
    text(fits[locate, ], labels=fits$brks[locate], cex = 0.7, offset = 0.2, pos = 2)

        #this is to extract the vector that identify brakpoints, for plotting in stratiplot
        zones <- fits[locate, ][,4]

    legend("bottomright", col = c(pcol, pcol, acol, "#3465a4"), pt.bg = c(pbg, pbg, abg, NA), pt.cex = pcex,
           pch = c(tpch, bpch, apch, NA), lty = c(NA, NA, NA, 1, 1), cex = 0.8,
           legend = c("Top","Bottom", "Breakpoint (cal yr BP)", "Total Phosphorous"), levels(a.factor(fits$traj)), bty = "n")

    dev.off()


## Make multittplot     
    par(mfrow=c(4,5))
    par(mar=c(2,2,2,2))
    
    lake <- c("Fondococha")
    lake_ages <- c("fondodocha")
    k=max(unique(sigClustFondo)) #number of groups in ccluster
    
    tt.lwd <- 1
    env_data <- env_data_red

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
  
  multittplot(ttracks[[lake]], lwd = tt.lwd, surf=NULL)
      
  coreCluster_lake <- coreCluster[[lake]]
  sigClust <- cutree(coreCluster_lake, k=k)
  
  agesClust <- as.numeric(coresList[[lake_ages]]$upper_age)
  
  fits <- as.data.frame(cbind(fitted(ttracks[[lake]], type = "passive"), 
                                  traj= sigClust, brks = agesClust))
      
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
  
  #calibrate () calls predict.gam
  newdat <- fits[,1:2]
  
  pred_space <- calibrate(surf, newdata = newdat, se.fit=TRUE)

  ordResult$modspace_fit[[i]] <- as.numeric(pred_space$fit)
  ordResult$modspace_fit_error[[i]] <- as.numeric(pred_space$se.fit)
  
  RMSE <- function(error) { sqrt(mean(error^2)) }
  ordResult$RMSE[[i]] <-  RMSE(surf$residuals)
  
  
}
    
  #name ordResult list  
  names(ordResult$fitted.values) <- colnames(env_data)
  names(ordResult$dev.expl) <- colnames(env_data)
  names(ordResult$s.pv) <- colnames(env_data)
  names(ordResult$grid) <- colnames(env_data)
  names(ordResult$modspace_fit) <- colnames(env_data)
  names(ordResult$modspace_fit_error) <- colnames(env_data)
  names(ordResult$RMSE) <- colnames(env_data)
  
  ordResult$dev.expl <- ordResult$dev.expl[order(ordResult$dev.expl,decreasing = TRUE)]
  ordResult$RMSE <- ordResult$RMSE[order(ordResult$RMSE,decreasing = TRUE)]
  
    # plot % deviance explained
    for (i in 1) {
      op <- par(mar = c(5,4,1,1) + 0.1)
      barplot(ordResult$dev.expl, cex.names=0.8, las=2, ylab="% deviance explained")
    }

    for (i in 1) {
      op <- par(mar = c(5,4,1,1) + 0.1)
      barplot(ordResult$RMSE, cex.names=0.8, las=2, ylab="RMSE")
    }
  
  #save list of results
  saveRDS(ordResult, file = "fondococha_gam_surface_space.rds")
  

###### Space-only model
    # extract predict values using PCA scrs from passive core ordination in ttracks
    lake <- "Llaviucu"
    lake <- "Yahuarcocha"
    lake <- "Pinan"
    
    sigClust <- sigClustLlav
    agesClust <- as.numeric(coresList$llaviucu$upper_age)
    
    sigClust <- sigClustYah
    agesClust <- as.numeric(coresList$yahuarcocha$upper_age)
    
    sigClust <- sigClustPin
    agesClust <- as.numeric(coresList$pinan$upper_age)
    
    # create sintetic data
    fits <- as.data.frame(cbind(fitted(ttracks[[lake]], type = "passive"), 
                    traj= sigClust, brks = agesClust))
    
    # using ordisurf()
    surf <- ordisurf(ttracks[[1]]$ord, env_data[,"P.season"],
                     method = "REML", scaling = ttracks[[1]]$scaling,
                     npoints = nrow(fits),
                     select = TRUE, plot = FALSE, knots = 20)
    
    newdat <- fits[,1:2]

    #calibrate () calls predict.gam
    pred_space <- calibrate(surf, newdata = newdat, se.fit=TRUE)
    fits$modspace_fit <- as.numeric(pred_space$fit)
    
    plot(fits$brks, fits$modspace_fit, col=fits$traj, type = "b")
    plot(fits$modspace_fit, fits$PC1, col=fits$traj, type="b")
    
  
   
  ##-- finer control using gam() and extracting scores
    scrs <- data.frame(scores(ttracks[[1]]$ord, display = "sites", scaling = 3))
    dat <- with(env_data, cbind(scrs, TP, K, Ca, Depth_avg, 
                                lake_depth_ratio, P.season, region_f=factor(env_data_lakes$region)))
    
    
    mod_spatial <- gam(TP ~ s(PC1, PC2, k = 10, bs="tp") + s(region_f, bs="re"),  
                       method="REML", data = dat, select=TRUE)
    
    sapply(mod_spatial$smooth, "[[", "label")
    summary(mod_spatial)
    
    
  ##-- newdata are PCA site scores from passive core ordination
    newdat <- fits[,1:2]
    newdat$region_f <- '0'
    
    # pass into 0s those covariates I don't have
    predict_temporal[c("region_f")] <- 0
    
    
  # extract predictions from model
    test <- predict(mod_spatial, newdata = newdat, type="terms", 
            terms = "s(PC1, PC2)", exclude = "s(region_f)")
    
    # mod_spatial <- gam(TP ~ s(PC1, PC2, k = 20, bs="tp") + s(Depth_avg) + ti(PC1, PC2, Depth_avg, d=c(2,1)), 
    #                    method="REML", data = dat, select=TRUE)
    

  ## tensor interactions allow multiple variables interact despite these having different scales
    tensor_mod <- gam(TP ~ te(PC1, PC2, P.season), 
                      data = dat, method = "REML", select=TRUE)
    
  # we can model only the interaction of two variables, and not their independent effects, 
    #which we estimate separately using ti()
    
    tensor_mod2 <- gam(TP ~ s(PC1, PC2) + s(lake_depth_ratio) + ti(PC1, PC2, lake_depth_ratio), 
                       data = dat, method = "REML", select=TRUE)
    plot(tensor_mod, pages=1)
    summary(tensor_mod)
    
    
########
# Proxy data
    llaviucu <- read.csv("llaviucu_diat_proxy_2RA.csv", row.names=1) %>%
      mutate(K_Ti = K/Ti)
    yahuarcocha <- read.csv("yahuarcocha_diat_proxy_2RA.csv", row.names = 1)
    umayo <- read.csv("umayo_diat_proxy_2RA.csv", row.names = 1)
    pinan <- read.csv("pinan_diat_proxy_2RA.csv", row.names = 1)
    fondococha <- read.csv("fondococha_diat_proxy_3RA.csv", row.names = 1)
    
    # Merge datasets
    proxyData <- analogue::join(llaviucu, yahuarcocha, pinan, umayo, fondococha, verbose = TRUE)
    

## Function to perform MULTIPLE PCA ordinations on diatom core data    
    doPCA <- function(i, cores, ...) {
      core <- cores[[i]]
      core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake"))] # drop year & depths vars
      core[is.na(core)] <- 0
      core <- core[, colSums(core) > 0]
      core <- decostand(core, method="hellinger") #Hellinger transform relative abundance data
      pca <- rda(core, scale=T, scaling=3)
      pca
    }

    #wrap up the function
    corePCA <- lapply(seq_along(coresList), doPCA, cores = coresList)
    names(corePCA) <- names(coresList)
    
    
    ## plot PCAs
    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2))
    
    for(i in seq_along(corePCA)) {
      op <- par(mar = c(5,4,1,1) + 0.1)
      plot(corePCA[[i]], type="n")
      #biplot(corePCA[[i]], scaling = 3)
      points(corePCA[[i]], col = "grey", pch = 19, cex = 0.7)
      title(names(corePCA[i]))
      par(op)
      
    }
    
    # plot PCA trajectories
    tt.lwd <- 1
    lwd <- 1
    
     
    for(i in seq_along(corePCA)) {
      PCAplot(corePCA[[i]], surf = NULL)
      title(names(corePCA[i]))
    }
    

    ## Plot species scores
    for (i in seq_along(corePCA)) {
      plot(summary(corePCA[[i]])$sites[,1], summary(corePCA[[i]])$sites[,2], cex=0.6, xlab="PC1", ylab="PC2", type="n", ylim=c(-1.6,1.4), xlim = c(-3, 3))
      abline(h=0, v=0, lty=2)
      points(summary(corePCA[[i]])$sites[,1], summary(corePCA[[i]])$sites[,2], cex=0.6, xlab="PC1", ylab="PC2", col="dark grey", pch=20 )
      
      
      loadPC1 <- order(abs(summary(corePCA[[i]])$species[,1]), decreasing=T)[1:8]
      loadPC2 <- order(abs(summary(corePCA[[i]])$species[,2]), decreasing=T)[1:8]
      
      specPC1 <- summary(corePCA[[i]])$species[,1]
      specPC2 <- summary(corePCA[[i]])$species[,2]
      spec <- c(loadPC1, loadPC2)
      arrows(x0=rep(0,24), y0=rep(0,24), x1=c(specPC1[spec]), y1=c(specPC2[spec]), col="dark red", length=0.07)
      text(x=c(specPC1[spec]), y=c(specPC2[spec]), labels=names(specPC1[spec]), cex=0.6)
      
      title(names(corePCA[i]))
    }
    

##-- Bin proxy data
proxyData$umayo <- NULL
    
binProxyData <- function(i, cores, ...) {
      core <- cores[[i]]
      #core <- core[ , which(colnames(core) %in% c("Fe_Mn", "Si_Ti", "K", "Ca", "Ti", "Si", "d13C", "C_N", "d18O"))] # drop year & depths vars
      #core <- core[, colSums(core) > 0]
      years <- round((proxyData[[i]]$upper_age*(-1)+1950), digits=0) #this is to get ages_AD
      
      core <- core[ , -which(colnames(core) %in% c("depth", "upper_age"))] # drop year & depths vars
      
      # Function for binning samples in each core (Alistair Seddon)
      source("/Volumes/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
      
      #make the age categories
      diatomBin1 <- binFunc(as.data.frame(core), as.numeric(years), 100, 0, 1700) ##
      diatomBin2 <- binFunc(as.data.frame(core), as.numeric(years), 30, 1700, 2000) ##
      
      #merge the two binning dataframes for North lakes
      rbind.data.frame(diatomBin1, diatomBin2)
}
    
#wrap up the function
binnedProxies <- lapply(seq_along(proxyData), binProxyData, cores = proxyData)
names(binnedProxies) <- names(proxyData)
   
save(binnedProxies, file="binnedProxies.RData")

#extract 
nams <- names(binnedProxies)
for (i in seq_along(binnedProxies)) {
  assign(paste0("binnedProxy_", nams[i]), binnedProxies[[i]])
}


# binnedproxies_df <- plyr::ldply(binnedProxies, data.frame)
# ages <- plyr::ldply(as.numeric(rownames(binnedProxies[[1]])), data.frame)

         
##-- Prepare proxy data
    #extract PCA scores
    scrs <- data.frame(scores(corePCA$llaviucu, display = "sites", scaling = 3))
    scrs <- data.frame(scores(corePCA$yahuarcocha, display = "sites", scaling = 3))
    scrs <- data.frame(scores(corePCA$pinan, display = "sites", scaling = 3))
    scrs <- data.frame(scores(corePCA$fondococha, display = "sites", scaling = 3))
    
    ## cbind with proxy data
    dat <- with(proxyData$llaviucu, cbind(scrs, Fe_Mn, Si_Ti, K_Ti))
      #transform data
      dat <- transform(dat, Si=log10(Si+0.25))
    
    dat <- with(proxyData$yahuarcocha, cbind(scrs, d13C, d18O, C_N))
      #transform data
      dat <- transform(dat, d13C=log10(d13C+10), d18O=sqrt(d18O), C_N=log10(C_N+0.25))
    
    dat <- with(proxyData$pinan, cbind(scrs, d13C, C_N))
    
    #fondodocha
    dat <- with(proxyData$fondococha, cbind(scrs, Fe_Mn, Si_Ti, K_Ti))
    

    # Plot correlations
    pairs(dat, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
    
    
    # drop PCA scrs 
    proxyVar <- dat[,-c(1,2)]
    
    ordResultProxy <- list()
    
    lake <- "pinan"
    lake_ages <- "Pinan"
    k=4
  
    for(i in 1:length(proxyVar)){
      
      proxy <- as.data.frame(
        cbind(dat, brks = coresList[[lake]]$upper_age,
              traj =  sigClust <- cutree(coreCluster[[lake_ages]], k=k)))
      proxy <- transform(proxy, elapsedTime=c(NA, diff(brks)))
      proxy <- proxy[-1,]
      
    
      surf <- ordisurf(corePCA[[lake]], proxyVar[,i], method = "REML", scaling=3,
                       select = TRUE, plot = TRUE, main=colnames(proxyVar[i]))
      ordResultProxy$dev.expl[[i]] <- summary(surf)$dev.expl
      ordResultProxy$fitted.values[[i]] <- fitted(surf)
      ordResultProxy$grid[[i]] <- surf$grid[1:3] 
      ordResultProxy$s.pv[[i]] <- summary(surf)[["s.pv"]]
      ordResultProxy$ages <- coresList[[lake]]$upper_age
      
      # model with gam()
      mod_temporal <- gam(proxyVar[,i][-1] ~ s(PC1, PC2, k = 10, bs="tp"), 
                          weights = elapsedTime / mean(elapsedTime),  
                          method="REML", data = proxy, select=TRUE)
      
      newdat <- proxy[,1:2]
      
      RMSE <- function(error) { sqrt(mean(error^2)) }
      ordResultProxy$RMSE[[i]] <-  RMSE(mod_temporal$residuals)
      
      #predict.gam
      pred_temporal <- predict(mod_temporal, newdata = newdat, se.fit=TRUE)
      
      ordResultProxy$proxyData <- proxy
      ordResultProxy$modtemp_fit[[i]] <- as.numeric(pred_temporal$fit)
      ordResultProxy$modtemp_fit_se[[i]] <- as.numeric(pred_temporal$se.fit)
      ordResultProxy$ages <- coresList[[lake]]$upper_age
}
    
    names(ordResultProxy$fitted.values) <- colnames(proxyVar)
    names(ordResultProxy$dev.expl) <- colnames(proxyVar)
    names(ordResultProxy$s.pv) <- colnames(proxyVar)
    names(ordResultProxy$grid) <- colnames(proxyVar)
    names(ordResultProxy$modtemp_fit) <- colnames(proxyVar)
    names(ordResultProxy$modtemp_fit_se) <- colnames(proxyVar)
    names(ordResultProxy$RMSE) <- colnames(proxyVar)
    
    ordResultProxy$dev.expl <- ordResultProxy$dev.expl[order(ordResultProxy$dev.expl,decreasing = TRUE)]
    
    # plot % deviance explained
    for (i in 1) {
      op <- par(mar = c(5,4,1,1) + 0.1)
      barplot(ordResultProxy$dev.expl, cex.names=0.8, las=2, ylab="% deviance explained")
    }
    
    #save list of results
    saveRDS(ordResultProxy, file = "fondococha_gam_surface_time.rds")
    
    
##### Temporal-only model
    # Llaviucu
    proxy <- as.data.frame(
      cbind(dat, brks = coresList$llaviucu$upper_age,
            traj = sigClustLlav))
    proxy <- transform(proxy, elapsedTime=c(NA, diff(brks)))
    proxy$zones <- zones
    proxy <- proxy[-1,]
    
    #Yahuarcocha
    proxy <- as.data.frame(
      cbind(dat, brks = coresList$yahuarcocha$upper_age,
            traj = sigClustYah))
    proxy <- transform(proxy, elapsedTime=c(NA, diff(brks)))
    proxy$zones <- zones
    proxy <- proxy[-1,]


    #Pinan
    proxy <- as.data.frame(
      cbind(dat, brks = coresList$pinan$upper_age,
            traj = sigClustPin))
    proxy <- transform(proxy, elapsedTime=c(NA, diff(brks)))
    proxy$zones <- zones
    proxy <- proxy[-1,]
    
    
    #Fondococha
    proxy <- as.data.frame(
      cbind(dat, brks = coresList$fondococha$upper_age,
            traj = sigClustFondo))
    proxy <- transform(proxy, elapsedTime=c(NA, diff(brks)))
    
    proxy$zones <- zones
    proxy <- proxy[-1,]
    
    
    # model with gam()
    mod_temporal <- gam(d13C ~ s(PC1, PC2, k = 10, bs="tp"), 
                        weights = elapsedTime / mean(elapsedTime),  
                        method="REML", data = proxy, select=TRUE, na.action = 'na.omit')
    
    summary(mod_temporal)
    plot(mod_temporal, page=1, scheme = 2)
    
    
    # make predictions from PCA scores axis 1 and 2
    newdat <- proxy[,1:2]
    
    #predict.gam
    pred_temporal <- predict(mod_temporal, newdata = newdat, se.fit=TRUE)
    proxy$modtemp_fit <- as.numeric(pred_temporal$fit)
    
    ## ggplot 
    plt <- ggplot(proxy, aes(x=modtemp_fit, y=PC1, shape=factor(zones), color=brks)) +
      geom_point() + 
      geom_smooth(method=lm, alpha=0.2, color="grey") +
      scale_color_gradient(high = "orange")+
      theme_bw() +
      xlab(expression(paste (delta^13, "C \u2030")))+
      #xlab("d13C")+
      ggtitle("Piñan") +
      labs(color = "cal yr BP", shape="Diatom zones")
    plt  
    
    ggsave("pinan_phase_plot.png", plt, height = 8, width = 10)  
    
    
    # R plain quick plot
    plot(proxy$brks, proxy$modtemp_fit, col=proxy$zones, type = "b")
    plot(proxy$modtemp_fit, proxy$PC1, col=proxy$zones, type = "b")
    
    plot(proxy$modtemp_fit, proxy$PC2, col=proxy$zones, type = "b")
    plot(proxy$modtemp_fit, proxy$PC1, type = "b",
         col=seq_palette[proxy[,"zones"]], pch=19)
    
    
   
    
    #plot stratigraphically
    plot(proxy$modtemp_fit, proxy$brks, xlab="GAM fitted d18O", ylab="Age cal yr BP", 
         type="b", main = "Yahuarcocha",
         ylim = rev(range(proxy$brks)), col=seq_palette[proxy[,"zones"]], pch=19)
    legend("topleft",legend = c("Zone 1", "Zone 2", "Zone 3", "Zone 4"), col=c(seq_palette[1],seq_palette[2], 
                                                                               seq_palette[3], seq_palette[4]),inset=0.010,cex=0.9, pch=19)

    
proxy$lake <- "yahuarcocha"    

    
## Compare space and temporal GAM fitted values
ordResult <- readRDS("fondococha_gam_surface_space.rds") 
ordResultProxy <- readRDS("fondococha_gam_surface_time.rds")    

ordResult <- readRDS("yahuarcocha_gam_surface_space.rds") 
ordResultProxy <- readRDS("yahuarcocha_gam_surface_time.rds")    

     
   #extract dataframes from list
   # fits are the df from each core
   space_fit <- data.frame(sapply(ordResult$modspace_fit,c))
   space_fit_err <- data.frame(sapply(ordResult$modspace_fit_error, c))
   
   space_fit$ages <- ordResult$ages

   time_fit <- data.frame(sapply(ordResultProxy$modtemp_fit,c))
   time_fit$ages <- sapply(ordResultProxy$proxyData$brks,c)
   time_fit$commtype <- sapply(ordResultProxy$proxyData$traj,c)
      time_fit$zones <- zones[-1]
   
   time_fit$PC1 <- sapply(ordResultProxy$proxyData[,1],c)
   time_fit$PC2 <- sapply(ordResultProxy$proxyData[,2],c)
   
   tidyFits <- time_fit %>%
     gather(key=variable, value=fit, -commtype, -zones, -ages, -PC1, -PC2)
   
   plot(time_fit$K_Ti, time_fit$PC2, xlab="GAM fitted Ki_Ti", ylab="Diatom PC1", 
        main = "", col=seq_palette[time_fit[,"commtype"]], pch=19)
   
   
   pltSpaceTime<-ggplot(tidyFits, aes(x=zones, y=fit, fill=as.factor(zones))) +
     facet_wrap(variable~., scales = "free") +
     geom_violin(trim = FALSE) +
     geom_jitter(height = 0, width = 0.1)+
     stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
     geom_boxplot(width=0.1)+
     scale_fill_brewer(palette="Blues") + theme_classic() +
     labs(x="Diatom Community types", y = "") +
     theme(legend.position = "none")
   pltSpaceTime
   
 

   # join space and temporal GAM fits
   space_time_data <- cbind(space_fit[-1,], time_fit)
   space_time_data_err <- cbind(space_fit_err[-1,], time_fit_err)
   
   #write.csv(space_time_data, "yahuarcocha_gam_space_time.csv")
   
   #pairs correlations
   pairs(space_time_data, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
   
   
   # make violin plots
   scaleFactor <- max(space_time_data$TP) / max(space_time_data$K_Ti)
   
   ggplot(data=space_time_data)+
     geom_violin(aes(x=zones, y=K_Ti, fill=as.factor(zones)))+
     geom_violin(aes(x=zones, y=TP*scaleFactor, fill=as.factor(zones)))+
     scale_y_continuous(name="K_Ti", sec.axis=sec_axis(~./scaleFactor, name="TP"))
     

 
   tidyFits <- space_time_data %>%
     gather(key=variable, value=fit, -commtype, -zones) %>%
     filter(variable %in% c("TP", "pH", "K_Ti"))
   
   ggplot(data=tidyFits, aes(x=zones, y=fit, fill=as.factor(variable)) +
            geom_violin() +
            geom_jitter(width=0.1,alpha=0.2)+
            xlab("")+ 
            facet_wrap(~variable))

   
   pltSpaceTime<-ggplot(tidyFits, aes(x=zones, y=fit, fill=as.factor(zones))) +
     facet_wrap(variable~., scales = "free") +
     geom_violin(trim = FALSE) +
     stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
     geom_boxplot(width=0.1)+
     geom_jitter(width=0.1,alpha=0.3)+
     scale_fill_brewer(palette="Blues") + theme_classic() +
     labs(x="Diatom Community types", y = "Fitted GAM") +
     theme(legend.position = "none")
   pltSpaceTime
   
   
   # plot
   library(GGally)
   
   cores_gam_space_time2 <- cores_gam_space_time %>% select(-age) %>% mutate(id_f=factor(cores_gam_space_time$id))
   
   ggpairs(cores_gam_space_time2, columns = c("TP", "Cond", "Si_Ti"), ggplot2::aes(colour=id_f),
           title = "",
           theme_bw()) 
   
   
   ggplot <- function(...) 
   ggplot2::ggplot(...) + scale_color_brewer(palette="Set1") + scale_fill_brewer(palette="Set1")
   unlockBinding("ggplot",parent.env(asNamespace("GGally")))
   assign("ggplot",ggplot,parent.env(asNamespace("GGally")))
   
   pltSpaceTime <- ggpairs(cores_gam_space_time2, columns = c("TP", "Cond"), 
           title = "", 
           upper = list(continuous = wrap("cor", size = 3)),
           lower = list(continuous = wrap("smooth",
                                          alpha = 0.3,
                                          size = 1)),
           mapping = aes(color = id_f) +
     theme_bw() +
     #stat_ellipse(geom = "polygon", alpha = 1/2, aes(fill = commtype)) +
   labs(title = "Llaviucu"))
   pltSpaceTime
   
   
   ggsave("Llaviucu_spacetime.png", pltSpaceTime, height = 8, width = 10)
   
  
   # Plot stratigraphically 
   #reverse temporal axis
   par(mfrow=c(1,2), mar=c(8,8,2,2))
   
   #plot
   library(viridis)
   seq_palette <- viridis(4)
   
   
   locate <- cumsum(rle(as.numeric(space_time_data$commtype))$lengths)
   #this is to extract the vector that identify brakpoints, for plotting in stratiplot
   zones <- space_time_data[locate, ][,"ages"]
   
   
   plot(space_time_data$TP, space_time_data$ages, xlab="GAM fitted TP", ylab="Age cal yr BP", 
        type="b", main = "Space",
        ylim = rev(range(space_time_data$ages)), col=seq_palette[space_time_data[,"commtype"]], pch=19)
   #for(i in 1:4) abline(h=zones[i], col=i+1, lty=2, lwd=2)
   arrows(x-sdev, y, x+sdev, y, length=0.05, angle=90, code=3)
   
   plotCI(space_time_data$TP,space_time_data$ages, uiw = space_time_data_err$TP,add=T,pch=NA,sfrac=10, col="blue")
   
   
   plot(space_time_data$Si_Ti, space_time_data$ages, xlab="GAM fitted Si_Ti", ylab="", type="b", 
        main = "Time",
        ylim = rev(range(space_time_data$ages)), col=seq_palette[space_time_data[,"commtype"]], pch=19)
   #for(i in 1:4) abline(h=zones[i], col=i+1, lty=2, lwd=2)
   
   
    
    

###
   llaviucu_space_time_data  <- read.csv("llaviucu_gam_space_time.csv", row.names = 1)
   yahuarcocha_space_time_data <- read.csv("yahuarcocha_gam_space_time.csv", row.names=1)
   pinan_space_time_data <- read.csv("pinan_gam_space_time.csv", row.names=1)
   fondococha_space_time_data <- read.csv("fondococha_gam_space_time.csv", row.names=1)
   
   #llaviucu
   ages <- llaviucu_space_time_data[,names(llaviucu_space_time_data) %in% c("ages")]
   llaviucu_space_time_data <- llaviucu_space_time_data[,!names(llaviucu_space_time_data) %in% c("commtype")]
   llaviucu_space_time_data <- llaviucu_space_time_data[,!names(llaviucu_space_time_data) %in% c("ages")]
   
   varBin <- binFunc(as.data.frame(llaviucu_space_time_data), ages, 30, -60, 1680) 
   varBin$age <- rownames(varBin)
   
   varBinInter <- na.approx(varBin, na.rm = TRUE) #do interpolation between adjacent samples
   llaviucu_gam_df_space <- as.data.frame(varBinInter)
   
   row.has.na <- apply(llaviucu_gam_df_space, 1, function(x){any(is.na(x))})
   sum(row.has.na)
   llaviucu_gam_df_space <- llaviucu_gam_df_space[!row.has.na,]
   
   llaviucu_gam_df_space$id <- "llaviucu"
   
   write.csv(llaviucu_gam_df_space, "llaviucu_gam_df_space_time_binned.csv")
   
   #yahuarcocha
   ages <- yahuarcocha_space_time_data[,names(yahuarcocha_space_time_data) %in% c("ages")]
   yahuarcocha_space_time_data <- yahuarcocha_space_time_data[,!names(yahuarcocha_space_time_data) %in% c("commtype")]
   yahuarcocha_space_time_data <- yahuarcocha_space_time_data[,!names(yahuarcocha_space_time_data) %in% c("ages")]
   
   varBin <- binFunc(as.data.frame(yahuarcocha_space_time_data), ages, 30, -60, 1680) 
   varBin$age <- rownames(varBin)
   
   varBinInter <- na.approx(varBin, na.rm = TRUE) #do interpolation between adjacent samples
   yahuarcocha_gam_df_space <- as.data.frame(varBinInter)
   
   row.has.na <- apply(yahuarcocha_gam_df_space, 1, function(x){any(is.na(x))})
   sum(row.has.na)
   yahuarcocha_gam_df_space <- yahuarcocha_gam_df_space[!row.has.na,]
   
   yahuarcocha_gam_df_space$id <- "yahuarcocha"
   
   write.csv(yahuarcocha_gam_df_space, "yahuarcocha_gam_df_space_time_binned.csv")
   
   #pinan
   ages <- pinan_space_time_data[,names(pinan_space_time_data) %in% c("ages")]
   pinan_space_time_data <- pinan_space_time_data[,!names(pinan_space_time_data) %in% c("commtype")]
   pinan_space_time_data <- pinan_space_time_data[,!names(pinan_space_time_data) %in% c("ages")]
   
   varBin <- binFunc(as.data.frame(pinan_space_time_data), ages, 30, -60, 1200) 
   varBin$age <- rownames(varBin)
   
   varBinInter <- na.approx(varBin, na.rm = TRUE) #do interpolation between adjacent samples
   pinan_gam_df_space <- as.data.frame(varBinInter)
   
   row.has.na <- apply(pinan_gam_df_space, 1, function(x){any(is.na(x))})
   sum(row.has.na)
   pinan_gam_df_space <- pinan_gam_df_space[!row.has.na,]
   
   pinan_gam_df_space$id <- "pinan"
   
   write.csv(pinan_gam_df_space, "pinan_gam_df_space_time_binned.csv")
   
   #fondococha
   ages <- fondococha_space_time_data[,names(fondococha_space_time_data) %in% c("ages")]
   fondococha_space_time_data <- fondococha_space_time_data[,!names(fondococha_space_time_data) %in% c("commtype")]
   fondococha_space_time_data <- fondococha_space_time_data[,!names(fondococha_space_time_data) %in% c("ages")]
   
   varBin <- binFunc(as.data.frame(fondococha_space_time_data), ages, 30, -60, 870) 
   varBin$age <- rownames(varBin)
   
   varBinInter <- na.approx(varBin, na.rm = TRUE) #do interpolation between adjacent samples
   fondococha_gam_df_space <- as.data.frame(varBinInter)
   fondococha_gam_df_space$id <- "fondococha"
   
   write.csv(fondococha_gam_df_space, "fondococha_gam_df_space_time_binned.csv")
   
   

   #bind rows
   composite <- bind_rows(llaviucu_gam_df_space, yahuarcocha_gam_df_space, 
                          pinan_gam_df_space, fondococha_gam_df_space)
   
    
   write.csv(composite, "cores_gam_space_time.csv")
   
    

###################################################################
## time series of BC dissimilarity
###################################################################

    #read data
    ts_BCdismilarity_north_lakes_bins <- read.csv("ts_BCdismilarity_north_lakes_bins.csv", row.names = 1)
    ts_BCdismilarity_central_lakes_bins <- read.csv("ts_BCdismilarity_central_lakes_bins.csv", row.names = 1)

    tsref_BCdismilarity_north_lakes_bins <- read.csv("tsref_BCdismilarity_north_lakes_bins.csv", row.names = 1)
    tsref_BCdismilarity_central_lakes_bins <- read.csv("tsref_BCdismilarity_central_lakes_bins.csv", row.names = 1)

    #this is with core binned data
    fossildata <- coreBinned_north$Llaviucu
    fossildata <- coreBinned_north$Yahuarcocha
    fossildata <- coreBinned_north$Pinan
    fossildata <- coreBinned_north$Fondococha

    fossildata <- coreBinned_central$Titicaca
    fossildata <- coreBinned_central$Umayo
    fossildata <- coreBinned_central$Lagunillas

      #this is with non-binned core  data
      fossildata <- north_diat$Llaviucu
      fossildata <- north_diat$Yahuarcocha
      fossildata <- north_diat$Pinan
      fossildata <- north_diat$Fondococha


    #Calculate species turnover
    A <- analogue::distance(fossildata, method="bray")
    
    A # How similar are different time steps?

    # create an indicator for all diagonals in the matrix
    d <- row(A) - col(A)
    d # We create an index identifying main diagonal and minor diagonals > we want diagonal 1 (1 position below main)

    # use split to group on these values
    diagonal<-split(A, d)
    diagonal # list, identifies vectors for each diagonal

    timeseriesBC<-unlist(diagonal["1"]) # select relevant one (diag = 1)
    #timeseriesBC<-unlist(diagonal["1"])/dif # select relevant one divided by age intervals (diag = 1)

    timeseriesBC # OK
    length(timeseriesBC)

    timeseriesBC_fromreference<-A[,ncol(A)] #take the oldest sample as reference to see the direction of change
    #timeseriesBC_fromreference<-A[,ncol(A)]/dif #take the oldest sample as reference divided by age intervals

    # plot time series
    par(mfrow=c(2,1))
    plot.ts(timeseriesBC) # time series of turnover rates
    plot.ts(timeseriesBC_fromreference) # time series of distance from the oldest sample


    #####
    # to plot
    Distance_core <- as.data.frame(timeseriesBC) #this is the ts of turnover rates
    Distance_fromreference <- as.data.frame(timeseriesBC_fromreference) #this is the ts from reference

      #LLAVIUCU
      # this is to add variable age and sigClust to the df
      Distance_core$age <- as.numeric(row.names(coreBinned_north$Llaviucu[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_north$Llaviucu))
      Distance_core$clust <- sigClustLlav[-1]
      
          #this is to add from non-binned data
          Distance_core$age <- as.numeric(north$llaviucu$upper_age[-1])
          Distance_fromreference$age <- as.numeric(north$llaviucu$upper_age)
          Distance_core$clust <- sigClustLlav[-1]
          Distance_fromreference$clust <- sigClustLlav
          
          # rename dataframes
          ts_llaviucu <- Distance_core
          tsref_llaviucu <- Distance_fromreference
          ts_llaviucu$id <- "llaviucu"
          tsref_llaviucu$id <- "llaviucu"
          
          
          test <- Distance_core %>% 
            group_by(clust) %>% 
            summarise(bd=mean(timeseriesBC))

          #Anderson dispersion method
          fossildata <- north_diat$Llaviucu
          age <- as.numeric(north$llaviucu$upper_age)
          intervals <- diff(north$llaviucu$upper_age)
          
            dist <- as.matrix(vegdist(fossildata, method="bray"))
            dist2 <- dist[row(dist) == col(dist) + 1] ## extract off-diagonal
            
            roc2 <- dist2 / intervals
            mod <- aov(roc2 ~ as.factor(sigClustLlav[-1]))
            summary(mod)
            (mod.HSD <- TukeyHSD(mod))
            plot(mod.HSD)
            
            
          dist <- vegdist(fossildata, method="bray")
          mod <- betadisper(dist, sigClustLlav, type="centroid")
          anova(mod)
          plot(mod)
          mod.HSD <- TukeyHSD(mod)
          plot(mod.HSD)
          
            
            ## plot roc2 against age
            plot(roc2 ~ head(age, -1), type = "l", ylab = "Rate of Change", xlab = "Age")
            
        
        
      #YAHUARCOCHA      
      Distance_core$age <- as.numeric(row.names(coreBinned_north$Yahuarcocha[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_north$Yahuarcocha))
          
          #this is to add from non-binned data
          Distance_core$age <- as.numeric(north$yahuarcocha$upper_age[-1])
          Distance_fromreference$age <- as.numeric(north$yahuarcocha$upper_age)
          Distance_core$clust <- sigClustYah[-1]
          Distance_fromreference$clust <- sigClustYah
          
          ts_yahuarcocha <- Distance_core
          tsref_yahuarcocha <- Distance_fromreference
          ts_yahuarcocha$id <- "yahuarcocha"
          tsref_yahuarcocha$id <- "yahuarcocha"
          
          test <- Distance_core %>% 
            group_by(clust) %>% 
            summarise(bd=sd(timeseriesBC))
          
          #Anderson dispersion method
          fossildata <- north_diat$Yahuarcocha
          age <- as.numeric(north$yahuarcocha$upper_age)
          intervals <- diff(north$yahuarcocha$upper_age)
          
            dist <- as.matrix(vegdist(fossildata, method="bray"))
            dist2 <- dist[row(dist) == col(dist) + 1] ## extract off-diagonal
            
            roc2 <- dist2 / intervals
            mod <- aov(roc2 ~ as.factor(sigClustYah[-1]))
            summary(mod)
            (mod.HSD <- TukeyHSD(mod))
            plot(mod.HSD)
            
          dist <- vegdist(fossildata, method="bray")
          mod <- betadisper(dist, sigClustYah, type = "centroid")
          anova(mod)
          plot(mod)
          (mod.HSD <- TukeyHSD(mod))
          plot(mod.HSD)
          
          ## plot roc2 against age
          plot(roc2 ~ head(age, -1), type = "l", ylab = "Rate of Change", xlab = "Age")
          

      #PINAN
      Distance_core$age <- as.numeric(row.names(coreBinned_north$Pinan[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_north$Pinan))
          
          #this is to add from non-binned data
          Distance_core$age <- as.numeric(north$pinan$upper_age[-1])
          Distance_fromreference$age <- as.numeric(north$pinan$upper_age)
          Distance_core$clust <- sigClustPin[-1]
          Distance_fromreference$clust <- sigClustPin
          
          ts_pinan <- Distance_core
          tsref_pinan <- Distance_fromreference
          ts_pinan$id <- "pinan"
          tsref_pinan$id <- "pinan"
          
          test <- Distance_core %>% 
            group_by(clust) %>% 
            summarise(bd=sd(timeseriesBC))
          
          #Anderson dispersion method
          #drop group 4 in Pinan because tephra
          sigClustPin_red<-sigClustPin[-(50:55)]
          
          fossildata <- north_diat$Pinan
          age <- as.numeric(north$pinan$upper_age)
          intervals <- diff(north$pinan$upper_age)
          
              dist <- as.matrix(vegdist(fossildata, method="bray"))
              dist2 <- dist[row(dist) == col(dist) + 1] ## extract off-diagonal
              
              roc2 <- dist2 / intervals
              roc2inf <- c(37,38,39,40,50,51) #infinite values
              roc2 <- roc2[is.finite(roc2)] 
              sigClustPin_red <- as.numeric(sigClustPin)[-roc2inf]
              
              mod <- aov(roc2 ~ as.factor(sigClustPin_red[-1]))
              summary(mod)
              (mod.HSD <- TukeyHSD(mod))
              plot(mod.HSD)
              
              ## plot roc2 against age
              plot(roc2 ~ head(age, -1), type = "l", ylab = "Rate of Change", xlab = "Age")
              
          
          dist <- vegdist(fossildata[-(50:55),], method="bray")

          mod <- betadisper(dist, sigClustPin_red, type = "centroid")
          anova(mod)
          plot(mod)
          (mod.HSD <- TukeyHSD(mod))
          plot(mod.HSD)
          
      #FONDOCOCHA    
      Distance_core$age <- as.numeric(row.names(coreBinned_north$Fondococha[- 1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_north$Fondococha))
          
          #this is to add from non-binned data
          Distance_core$age <- as.numeric(north$fondococha$upper_age[-1])
          Distance_fromreference$age <- as.numeric(north$fondococha$upper_age)
          Distance_core$clust <- sigClustFondo[-1]
          Distance_fromreference$clust <- sigClustFondo
          
          
          ts_fondococha <- Distance_core
          tsref_fondococha <- Distance_fromreference
          ts_fondococha$id <- "fondococha"
          tsref_fondococha$id <- "fondococha"
          
          
          test <- Distance_core %>% 
            group_by(clust) %>% 
            summarise(bd=sd(timeseriesBC))
          
          #Anderson dispersion method
          fossildata <- north_diat$Fondococha
          age <- as.numeric(north$fondococha$upper_age)
          intervals <- diff(north$fondococha$upper_age)
          
            dist <- as.matrix(vegdist(fossildata, method="bray"))
            dist2 <- dist[row(dist) == col(dist) + 1] ## extract off-diagonal
            
            roc2 <- dist2 / intervals
            roc2inf <- c(82,83) #infinite values
            roc2 <- roc2[is.finite(roc2)] 
            sigClustFondo_red <- as.numeric(sigClustFondo)[-roc2inf]

            mod <- aov(roc2 ~ as.factor(sigClustFondo_red[-1]))
            summary(mod)
            (mod.HSD <- TukeyHSD(mod))
            plot(mod.HSD)
            
          
          dist <- vegdist(fossildata, method="bray")
          mod <- betadisper(dist, sigClustFondo, type = "centroid")
          anova(mod)
          plot(mod)
          (mod.HSD <- TukeyHSD(mod))
          plot(mod.HSD)
          
      ts_cores_df_north <- rbind(ts_llaviucu,
                                     ts_yahuarcocha,
                                     ts_pinan,
                                     ts_fondococha)
      write.csv(ts_cores_df_north, "ts_BCdismilarity_north_lakes.csv")
          
      tsref_cores_df_north <- rbind(tsref_llaviucu,
                                    tsref_yahuarcocha,
                                    tsref_pinan,
                                    tsref_fondococha)
      write.csv(tsref_cores_df_north, "tsref_BCdismilarity_north_lakes.csv")
      
              
      # this is to add variable age to the df
      Distance_core$age <- as.numeric(row.names(coreBinned_central$Titicaca[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_central$Titicaca))

      # this is to add variable age to the df
      Distance_core$age <- as.numeric(row.names(coreBinned_central$Umayo[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_central$Umayo))

      # this is to add variable age to the df
      Distance_core$age <- as.numeric(row.names(coreBinned_central$Lagunillas[-1,]))
      Distance_fromreference$age <- as.numeric(row.names(coreBinned_central$Lagunillas))

        #something odd with Lagunillas
        Distance_core$age <- as.numeric(row.names(lagunillasBin[-c(1,2),]))
        Distance_fromreference$age <- as.numeric(row.names(lagunillasBin[-1,]))


   
    ts_titicaca <- Distance_core
    tsref_titicaca <- Distance_fromreference

    ts_umayo <- Distance_core
    tsref_umayo <- Distance_fromreference

    ts_lagunillas <- Distance_core
    tsref_lagunillas <- Distance_fromreference


    #assign id to each ts of BC
    ts_titicaca$id <- "titicaca"
    ts_umayo$id <- "umayo"
    ts_lagunillas$id <- "lagunillas"

    #assign id to each ts from reference of BC
    tsref_titicaca$id <- "titicaca"
    tsref_umayo$id <- "umayo"
    tsref_lagunillas$id <- "lagunillas"

    ts_cores_df_north <- rbind(ts_llaviucu,
                               ts_yahuarcocha,
                               ts_pinan,
                               ts_fondococha)
    write.csv(ts_cores_df_north, "ts_BCdismilarity_north_lakes_bins.csv")

    ts_cores_df_central <- rbind(ts_titicaca,
                                 ts_umayo,
                                 ts_lagunillas)
    write.csv(ts_cores_df_central, "ts_BCdismilarity_central_lakes_bins.csv")

    tsref_cores_df_north <- rbind(tsref_llaviucu,
                                  tsref_yahuarcocha,
                                  tsref_pinan,
                                  tsref_fondococha)
    write.csv(tsref_cores_df_north, "tsref_BCdismilarity_north_lakes_bins.csv")

    tsref_cores_df_central <- rbind(tsref_titicaca,
                                    tsref_umayo,
                                    tsref_lagunillas)
    write.csv(tsref_cores_df_central, "tsref_BCdismilarity_central_lakes_bins.csv")


    BC_plot <- ggplot(ts_BCdismilarity_north_lakes_bins, aes(x=age, y=timeseriesBC, color=id))+
      geom_line()+
      #facet_grid(id~., scales = "free") +
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("BC dissimilarity") +
      theme_bw()

    #ggsave("BC_NorthAndes.png", BC_plot, height = 8, width = 10)

    BCref_plot <- ggplot(tsref_BCdismilarity_north_lakes_bins, aes(x=age, y=timeseriesBC_fromreference,
                                                                   color=id))+
      #facet_grid(id~., scales = "free") +
      geom_line()+
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("BC dissimilarity from reference") +
      theme_bw()

    #ggsave("BCref_NorthAndes.png", BCref_plot, height = 8, width = 10)


    library(cowplot)

    plt <- plot_grid(BC_plot, BCref_plot, align="v", axis="tb", nrow = 2,
                     rel_widths = c(1,2))
    ggsave("BC_NorthAndes_bins.png", plt, height = 8, width = 10)

    #plot_grid(ts_diss, ts_ref_diss, align="v", axis="tb", nrow = 2, rel_widths = c(1,2))
    
    

    ###################################################################
    ## Decompose BC time-series in replacement and richness components
    ###################################################################

    library(adespatial)

    A <- beta.div.comp(fossildata, coef = "S", quant = TRUE)
    A # How similar are different time steps?


    Replmatrix <- as.matrix(A$repl)
    Nesmatrix <- as.matrix(A$rich)

    # create an indicator for all diagonals in the matrix
    d <- row(Replmatrix) - col(Replmatrix)
    d <- row(Nesmatrix) - col(Nesmatrix)

    d # We create an index identifying main diagonal and minor diagonals > we want diagonal 1 (1 position below main)

    # use split to group on these values
    diagonal<-split(Replmatrix, d)
    diagonal<-split(Nesmatrix, d)


    tsBC_repl<-unlist(diagonal["1"]) # select relevant one (diag = 1)
    tsBC_nedt<-unlist(diagonal["1"]) # select relevant one (diag = 1)

    #Replmatrix_fromreference<-Replmatrix[,ncol(Replmatrix)] #take the oldest sample as reference to see the direction of change
    #Nesmatrix_fromreference<-Nesmatrix[,ncol(Nesmatrix)] #take the oldest sample as reference to see the direction of change


    # plot time series
    par(mfrow=c(2,1))
    plot.ts(tsBC_repl) # time series of replacement
    plot.ts(tsBC_nedt) # time series of nestedness


    #####
    # to plot
    # Replacement
    Distance_core_repl <- as.data.frame(tsBC_repl) #this is the ts of turnover rates

    Distance_core_repl$age <- as.numeric(row.names(coreBinned_north$Llaviucu[-1,]))
    Distance_core_repl$age <- as.numeric(row.names(coreBinned_north$Yahuarcocha[-1,]))
    Distance_core_repl$age <- as.numeric(row.names(coreBinned_north$Pinan[-1,]))
    Distance_core_repl$age <- as.numeric(row.names(coreBinned_north$Fondococha[-1,]))

    Distance_core_repl$age <- as.numeric(row.names(coreBinned_central$Titicaca[-1,]))
    Distance_core_repl$age <- as.numeric(row.names(coreBinned_central$Umayo[-1,]))
    Distance_core_repl$age <- as.numeric(row.names(coreBinned_central$Lagunillas[-1,]))


    #Nestedness
    Distance_core_nest <- as.data.frame(tsBC_nedt) #this is the ts of nestedness
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_north$Llaviucu[-1,]))
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_north$Yahuarcocha[-1,]))
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_north$Pinan[-1,]))
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_north$Fondococha[-1,]))

    Distance_core_nest$age <- as.numeric(row.names(coreBinned_central$Titicaca[-1,]))
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_central$Umayo[-1,]))
    Distance_core_nest$age <- as.numeric(row.names(coreBinned_central$Lagunillas[-1,]))

###
    # rename dataframes
    ts_repl_llaviucu <- Distance_core_repl
    ts_nest_llaviucu <- Distance_core_nest

    ts_repl_yahuarcocha <- Distance_core_repl
    ts_nest_yahuarcocha <- Distance_core_nest

    ts_repl_pinan <- Distance_core_repl
    ts_nest_pinan <- Distance_core_nest

    ts_repl_fondococha <- Distance_core_repl
    ts_nest_fondococha <- Distance_core_nest

    ts_repl_titicaca <- Distance_core_repl
    ts_nest_titicaca <- Distance_core_nest

    ts_repl_umayo <- Distance_core_repl
    ts_nest_umayo <- Distance_core_nest

    ts_repl_lagunillas <- Distance_core_repl
    ts_nest_lagunillas <- Distance_core_nest

###
    #assign id to each ts of BC
    ts_repl_llaviucu$id <- "llaviucu"
    ts_nest_llaviucu$id <- "llaviucu"

    ts_repl_yahuarcocha$id <- "yahuarcocha"
    ts_nest_yahuarcocha$id <- "yahuarcocha"

    ts_repl_pinan$id <- "pinan"
    ts_nest_pinan$id <- "pinan"

    ts_repl_titicaca$id <- "titicaca"
    ts_nest_titicaca$id <- "titicaca"

    ts_repl_umayo$id <- "umayo"
    ts_nest_umayo$id <- "umayo"

    ts_repl_lagunillas$id <- "lagunillas"
    ts_nest_lagunillas$id <- "lagunillas"

####
    ts_repl_df_north <- rbind(ts_repl_llaviucu,
                              ts_repl_yahuarcocha,
                              ts_repl_pinan)


    ts_nestd_df_north <- rbind(ts_nest_llaviucu,
                               ts_nest_yahuarcocha,
                               ts_nest_pinan)


    ts_repl_df_central <- rbind(ts_repl_titicaca,
                                ts_repl_umayo,
                                ts_repl_lagunillas)
    write.csv(ts_repl_df_central, "ts_repl_central_lakes.csv")


    ts_nestd_df_central <- rbind(ts_nest_titicaca,
                                 ts_nest_umayo,
                                 ts_nest_lagunillas)
    write.csv(ts_nestd_df_central, "ts_nest_central_lakes.csv")

#####
        BC_repl_plt <- ggplot(ts_repl_df_central, aes(x=age, y=tsBC_repl, color=id))+
          geom_line()+
          theme(legend.position="none") +
          xlab("Cal yr BP") + ylab("BC species replacement") +
          theme_bw()

        BC_nest_plt <- ggplot(ts_nestd_df_central, aes(x=age, y=tsBC_nedt, color=id))+
          geom_line()+
          theme(legend.position="none") +
          xlab("Cal yr BP") + ylab("BC nestedness") +
          theme_bw()

        library(cowplot)

        plt <- plot_grid(BC_repl_plt, BC_nest_plt, align="v", axis="tb", nrow = 2, rel_widths = c(1,2))
        ggsave("BC_podani_centralAndes.png", plt, height = 8, width = 10)



    # this is to add variable age to the df
    Distance_core$age <- as.numeric(coresList$llaviucu$upper_age[-1])
    Distance_fromreference$age <- as.numeric(coresList$llaviucu$upper_age)

    Distance_core$age <- as.numeric(coresList$yahuarcocha$upper_age[-1])
    Distance_fromreference$age <- as.numeric(row.names(coresList$yahuarcocha))



###################################################################
## Distance-based core trajectories (Lamothe et al 2017 Ecosphere)
###################################################################

   
#core trajectories are the fits data.frame from timetrack analysis: PCA axis 1 and 2 scores
##This is to calculate euclidean distances from consecutive samples standardized by age interval

# what if I calculate velocity of change using PCA scrs from core ordinations and not from timetrack?  
# same results than using PCA passive scores
    
    # scrs$brks <- fits$brks
    # scrs$traj <- fits$traj
    # 
    #     fits <- scrs
        
### 
   time <- fits$brks
    
   time <- ordResult$fits$brks    
   fits <- ordResult$fits  
        

        #this is to calculate time lag between observations
        dif <- c(NA, diff(time)) #NA

        time <- length(fits$brks)

        dist_matrix <- matrix(NA, nrow = nrow(fits), ncol = 4, byrow = FALSE, dimnames = NULL)
        colnames(dist_matrix) <- c("PC1", "PC2", "ElapsedTime", "Distance")
        dist_matrix[,1] <-fits$PC1
        dist_matrix[,2] <-fits$PC2
        dist_matrix[,3] <-c(NA, diff(fits$brks))

     
        for (x in 1:(length(dist_matrix))){
          for (i in 1:time){
            dist_matrix[i,4]<-sqrt((dist_matrix[i,1]-mean(dist_matrix[1:nrow(fits),1]))^2+
                                     (dist_matrix[i,2]-mean(dist_matrix[1:nrow(fits),2]))^2)/dist_matrix[i,3]
          }
        }
        
        Distance_core <- as.data.frame(dist_matrix)
        Distance_core$age <- fits$brks
        Distance_core$traj <- fits$traj
        
        plot(Distance_core$age, Distance_core$Distance, type = "l", col=as.factor(Distance_core$traj))
        
        
        # Reference lakes
        distances_fondococha <- Distance_core
        distances_fondococha$lake <- "fondococha"
        distances_pinan <- Distance_core
        distances_pinan$lake <- "pinan"
        reference_cumulative_dist <- rbind(distances_fondococha, distances_pinan)

        write.csv(reference_cumulative_dist, "reference_cumulative_dist.csv")
        
        
        
        
        # Impactes lakes
        distances_llaviucu <- Distance_core
        distances_llaviucu$lake <- "llaviucu"
        
        distances_yahuarcocha <- Distance_core
        distances_yahuarcocha$lake <- "yahuarcocha"
        
        impacted_cumulative_dist <- rbind(distances_llaviucu, distances_yahuarcocha)
        impacted_cumulative_dist_binned <- rbind(distances_llaviucu, distances_yahuarcocha)
        
        write.csv(impacted_cumulative_dist_binned, "impacted_cumulative_dist_binned.csv")
        
        
    # This is to calculate distances from historical baseline, assumed to be the oldest clustered samples
    # or first observation (nrows)
        
        time <- length(fits$brks)

        dist_matrix <- matrix(NA, nrow = nrow(fits), ncol = 4, byrow = FALSE, dimnames = NULL)
        colnames(dist_matrix) <- c("PC1", "PC2", "ElapsedTime", "Distance")
        dist_matrix[,1] <-fits$PC1
        dist_matrix[,2] <-fits$PC2
        dist_matrix[,3] <-c(NA, diff(fits$brks))

        for (x in 1:(length(dist_matrix))){
          for (i in 1:time){
            dist_matrix[i,4]<-sqrt((dist_matrix[i,1]-mean(dist_matrix[nrow(fits),1]))^2+
                                     (dist_matrix[i,2]-mean(dist_matrix[nrow(fits),2]))^2)/dist_matrix[i,3]
          }
        }


        Distance_core <- as.data.frame(dist_matrix)
        Distance_core$age <- fits$brks
        Distance_core$traj <- fits$traj

        plot(Distance_core$age, Distance_core$Distance, type = "l", col=as.factor(Distance_core$traj))

        
        # reference lakes
        distances_pinan_from_reference <- Distance_core
        distances_pinan_from_reference$lake <- "pinan"
        
        distances_fondococha_from_reference <- Distance_core
        distances_fondococha_from_reference$lake <- "fondococha"
        
        reference_baseline_dist<- rbind(distances_fondococha_from_reference, 
                                         distances_pinan_from_reference)
        write.csv(reference_baseline_dist, "reference_baseline_dist_binned.csv")
        
        # Impactes lakes
        distances_llaviucu_from_reference <- Distance_core
        distances_llaviucu_from_reference$lake <- "llaviucu"
        
        
        distances_yahuarcocha_from_reference <- Distance_core
        distances_yahuarcocha_from_reference$lake <- "yahuarcocha"
        
        impacted_reference_dist_binned <- rbind(distances_llaviucu_from_reference, 
                                           distances_yahuarcocha_from_reference)
        write.csv(impacted_reference_dist_binned, "impacted_baseline_dist_binned.csv")
        

        distances_cores_successive <- rbind(distances_pinan, distances_yahuarcocha, distances_fondococha, distances_llaviucu)
        write.csv(distances_cores_successive, "distances_cores_successive.csv")
        
        
        distances_cores_cumulative <- rbind(distances_pinan_from_reference, distances_yahuarcocha_from_reference, 
                                            distances_fondococha_from_reference, distances_llaviucu_from_reference)
        write.csv(distances_cores_cumulative, "distances_cores_cumulative.csv")
        
        
# Define baselines, boundary of distances traveled in multivariate space by all fossil diatom assemblages
        coreDiat <- new
        coreDiat_var <- coreDiat[ , which(names(coreDiat) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars

        coreDiat <- coreDiat[ , -which(names(coreDiat) %in% c("depth","upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars

        coreDiat[is.na(coreDiat)] <- 0

        # PCA cores
            coreDiat_hell <- decostand(coreDiat, method="hellinger") #Hellinger transform relative abundance data
            coreDiat_pca <- rda(coreDiat_hell, scale=T)

            scrs <- scores(coreDiat_pca, choices = 1:2, display = "sites")
            PCAscrs_cores <- as.data.frame(scrs)

            PCAscrs_cores$upper_age <- new$upper_age
            newdata <- PCAscrs_cores[order(PCAscrs_cores$upper_age),]

            PCAscrs_cores <- newdata %>%
              filter(!upper_age == 0)

            # This is to calculate distances from consecutive samples standardized by age interval
            time <- length(PCAscrs_cores$upper_age)

            dist_matrix <- matrix(NA, nrow = nrow(PCAscrs_cores), ncol = 4, byrow = FALSE, dimnames = NULL)
            colnames(dist_matrix) <- c("PC1", "PC2", "ElapsedTime", "Distance")
            dist_matrix[,1] <-PCAscrs_cores$PC1
            dist_matrix[,2] <-PCAscrs_cores$PC2
            dist_matrix[,3] <-c(NA, diff(PCAscrs_cores$upper_age))

            for (x in 1:(length(dist_matrix))){
              for (i in 1:time){
                dist_matrix[i,4]<-sqrt((dist_matrix[i,1]-mean(dist_matrix[1:nrow(fits),1]))^2+
                                         (dist_matrix[i,2]-mean(dist_matrix[1:nrow(fits),2]))^2)/dist_matrix[i,3]
              }
            }

            Distances.cores <- as.data.frame(dist_matrix)
            Distances.cores$age <- PCAscrs_cores$upper_age


            Distances.cores.summary<-summarySE(Distances.cores, measurevar="Distance", groupvars=c("age"))


            ggplot(data=Distances.cores.summary, aes(x=age, y=mean)) +
              geom_ribbon(aes(ymin=0,ymax=mean+(ci)), fill="grey70") +
              geom_point(colour="black", shape=20, size=2) +
              geom_line() +
              geom_path(data=Distance_core, aes(x=age, y=Distance)) +
              geom_point(data=Impacted.cumulative.dist, aes(x=Year, y=Distance, colour=Community)) +
              ylab("Distance to Baseline Centroid") +
              ylim(0,1)+
              theme_bw()+
              theme(axis.text = element_text(colour = "black", size=rel(1.5)),
                    axis.title.y = element_text(colour = "black", size=rel(1.5)),
                    axis.title.x = element_text(colour = "black", size=rel(1.5)),
                    axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    legend.position = "none")






##############################################
############### OLD stuff ####################


    ##calculate distances from each core to all modern samples
    #A matrix of dissimilarities where columns are the samples in training set and the rows the samples in the cores
    hellingerTrans <- function(i, cores, ...) {
      core <- cores[[i]]
      core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
      core <- decostand(core, method="hellinger")
      return(core)
    }

    cores_dist <- lapply(seq_along(df), hellingerTrans, cores=df)
    #name list elements
    names(cores_dist) <- c("Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

    # distances from each core to all modern samples
    dist_llaviucu <- distance(cores_dist$Llaviucu, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$llaviucu$upper_age))
    roc_llaviucu <- dist_llaviucu/dif
    roc_llaviucu <- as.data.frame(roc_llaviucu)
    roc_llaviucu$age <- as.numeric(coresList$llaviucu$upper_age)

    roc_llaviucu$dif <- dif

    roc_llaviucu <- roc_llaviucu %>%
      gather(key = sites, value = dist, -age) %>%
      group_by(age) %>%
      summarise(
        MEAN = mean(dist),
        SD = sd((dist), na.rm = TRUE),
        N = sum(!is.na(dist)),
        upper_limit = MEAN + SD/sqrt(N),
        lower_limit = MEAN - SD/sqrt(N)
      )

    distances_plot <- ggplot(roc_llaviucu, aes(x=age, y=MEAN, group=1))+
      geom_line(aes(x=age, y=MEAN, group=1), size=1)+
      geom_ribbon(aes(ymin=lower_limit,ymax=upper_limit),color="grey",alpha=0.4) +
      theme_light(base_size = 16) +
      #scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +
      #facet_grid(id_f~.) +
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("Distance to reference")


    #distances between consecutive samples
    dist_llaviucu <- distance(cores_dist$Llaviucu, method = "bray")
    dif <- c(NA, diff(coresList$llaviucu$upper_age))
    SCDcrop<-dist_llaviucu[-1,]
    scd<- diag(SCDcrop)
    roc <- scd/dif

    roc_llaviucu <- as.data.frame(roc)
    roc_llaviucu$age <- as.numeric(coresList$llaviucu$upper_age)

    distances_plot <- ggplot(roc_llaviucu, aes(x=age, y=roc, group=1))+
      geom_line(aes(x=age, y=roc, group=1), size=1)+
      #geom_ribbon(aes(ymin=lower_limit,ymax=upper_limit),color="grey",alpha=0.4) +
      theme_light(base_size = 16) +
      #scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +
      #facet_grid(id_f~.) +
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("Distance to reference")

    #Albert's code
    vector_data_disimi<- unmatrix(dist_llaviucu,byrow=F)
    lower_triangle<-lower.tri(dist_llaviucu)
    vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
    disimi_result<-vector_data_disimi[vector_data_triangle]





    dist_yahuarcocha <- distance(cores_dist$Yahuarcocha, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$yahuarcocha$upper_age))
    roc_yahuarcocha <- dist_yahuarcocha/dif
    roc_yahuarcocha <- as.data.frame(roc_yahuarcocha)
    roc_yahuarcocha$age <- as.numeric(coresList$yahuarcocha$upper_age)


    #distances between consecutive samples
    dist_yahuarcocha <- distance(cores_dist$Yahuarcocha, method = "bray")
    dif <- c(NA, diff(coresList$yahuarcocha$upper_age))
    SCDcrop<-dist_yahuarcocha[-1,]
    scd<- diag(SCDcrop)
    roc <- scd/dif

    roc_yahuarcocha <- as.data.frame(roc)
    roc_yahuarcocha$age <- as.numeric(coresList$yahuarcocha$upper_age)

    distances_plot <- ggplot(roc_yahuarcocha, aes(x=age, y=roc, group=1))+
      geom_line(aes(x=age, y=roc, group=1), size=1)+
      #geom_ribbon(aes(ymin=lower_limit,ymax=upper_limit),color="grey",alpha=0.4) +
      theme_light(base_size = 16) +
      #scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +
      #facet_grid(id_f~.) +
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("Distance")


    dist_titicaca <- distance(cores_dist$Titicaca, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$titicaca$upper_age))
    roc_titicaca <- dist_titicaca/dif
    roc_titicaca <- as.data.frame(roc_titicaca)
    roc_titicaca$age <- as.numeric(coresList$titicaca$upper_age)




    dist_umayo <- distance(cores_dist$Umayo, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$umayo$upper_age))
    roc_umayo <- dist_umayo/dif
    roc_umayo <- as.data.frame(roc_umayo)
    roc_umayo$age <- as.numeric(coresList$umayo$upper_age)



    dist_pinan <- distance(cores_dist$Pinan, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$pinan$upper_age))
    roc_pinan <- dist_pinan/dif
    roc_pinan <- as.data.frame(roc_pinan)
    roc_pinan$age <- as.numeric(coresList$pinan$upper_age)



    dist_triumfo <- distance(cores_dist$Triumfo, cores_dist$trainingset, method = "bray")
    dif <- c(NA, diff(coresList$triumfo$upper_age))
    roc_triumfo <- dist_triumfo/dif
    roc_triumfo <- as.data.frame(roc_triumfo)
    roc_triumfo$age <- as.numeric(coresList$triumfo$upper_age)



    #plot
    roc_llaviucu$id <- "llaviucu"
    roc_yahuarcocha$id <- "yahuarcocha"
    roc_titicaca$id <- "titicaca"
    roc_umayo$id <- "umayo"
    roc_pinan$id <- "pinan"
    roc_triumfo$id <- "triumfo"

    df_merged <- bind_rows(roc_llaviucu, roc_yahuarcocha, roc_titicaca, roc_umayo, roc_pinan, roc_triumfo)

    df_tidy <- gather(df_merged, sites, dist, -id,-age)
    df_tidy <- unite(df_tidy, unique_id, c(id, sites), sep="_", remove = FALSE)

    df_tidy_mean <- df_tidy %>%
      filter(!is.na(dist)) %>%
      group_by(age, id) %>%
      summarise(n = n(),
                mean = mean(dist),
                median = median(dist),
                sd = sd(dist)) %>%
      mutate(sem = sd / sqrt(n - 1),
             CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
             CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

    # this is to manually sort lakes
    df_tidy_mean$id_f = factor(df_tidy_mean$id, levels=c("triumfo", "pinan", "yahuarcocha", "llaviucu", "umayo", "titicaca"))

    color_list <- c("darkgoldenrod", "limegreen", "turquoise2", "blue", "red", "yellow")

    distances_plot <- ggplot(df_tidy_mean, aes(x=age, y=mean, color = id, group=1))+
      geom_line(aes(x=age, y=mean, color=id, group=1), size=1)+
      geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=id),color="grey",alpha=0.4) +
      theme_light(base_size = 16) +
      scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +
      facet_grid(id_f~., scales = "free") +
      geom_smooth(method = "loess") +
      theme(legend.position="none") +
      xlab("Cal yr BP") + ylab("Distance to reference")

    ggsave("Distances_to_reference_plot.png", distances_plot, height = 8, width = 10)


    llaviucu <- df_tidy_mean %>%
      filter(id == "triumfo")


    plt <- ggplot(llaviucu, aes(x=age, y=mean))+
      geom_line(aes(x=age, y=mean), size=1)+
      geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, group=1),color="grey",alpha=0.4) +
      theme_light(base_size = 16) +
      scale_x_reverse() +
      geom_smooth(method = "loess") +
      scale_fill_manual(values=color_list) + scale_color_manual(values=color_list) +
      theme(legend.position="none") +
      xlab("Age (cal yr BP)") + ylab("Distance to reference")

    ggsave("llaviucu_distance.png", plt, height = 8, width = 10)

    # fit a gam
    mod1 <- gam(mean ~ s(age), family = gaussian(link="identity"),
                data = llaviucu, method = "REML")


    all_distances <- all_distances %>%
      filter(lake == "llaviucu") %>%
      mutate(AgeCE = age*(-1)+1950) %>%
      filter(AgeCE >= 0)




    # BASELINE DISTANCES
    # join cores_merged and trainingset

    cores_merged_sorted_by_age <- cores_merged[order(cores_merged$upper_age),]
    cores_training_joined <- join(cores_merged_sorted_by_age, training)

    #drop variables to transform hellinger relaative abundances
    remove <- function(i, cores, ...) {
      core <- cores[[i]]
      core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
      core <- decostand(core, method="hellinger")
      return(core)
    }

    cores_training_clean <- lapply(seq_along(cores_training_joined), remove, cores=cores_training_joined)
    #name list elements
    names(cores_training_clean) <- c("cores", "trainingset")

    #calculate distances from each core sample to all samples in the training set
    dist <- distance(cores_training_clean$cores, cores_training_clean$training, method = "bray")

    #calculate distances between consecutive samples
    cores_dist <- cores_merged_sorted_by_age[, -which(names(cores_merged_sorted_by_age) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))]
    cores_dist[is.na(cores_dist)] <- 0
    cores_hell <- decostand(cores_dist, method="hellinger")
    #rownames(cores_hell) <- as.numeric(cores_merged_sorted_by_age$upper_age)

    dist <- distance(cores_hell, method = "bray")
    #dist_ord <- sort_dist_mat(dist)

    # unroll the dissimilarity matrix
    vector_data_disimi<- unmatrix(dist,byrow=F)
    lower_triangle<-lower.tri(dist)
    vector_data_triangle<-unmatrix(lower_triangle, byrow=F)
    disimi_result<-vector_data_disimi[vector_data_triangle]

    #standardized distances by age intervals between samples
    dif <- c(NA, diff(cores_merged_sorted_by_age$upper_age))

    rownames(dist) <- as.numeric(cores_merged_sorted_by_age$upper_age)

    # divide distances by age intervals between samples (not sure if this is the correct operation)
    dist_crop<-dist[-1,]
    dist_diag<- diag(dist_crop)

    dif <- c(NA, diff(cores_merged_sorted_by_age$upper_age))[-1]
    roc <- dist_diag/dif
    plot(roc)

    mean <- apply(scd, 1, mean)



