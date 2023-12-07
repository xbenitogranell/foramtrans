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
library(mgcv)

## Read in forams counts
forams <- read.csv("datasets/S2/S2_counts_v2.csv", sep=";")[-1] %>%
  mutate(depth=depth*100) #depth in cm
str(forams)

# Read in forams taxa names groups
changes <- read.csv("datasets/S2/nms_taxa_groups.csv")

# Read in age-depth model
ages <- read.table("Bacon_runs/S2_core_v2/S2_core_v2_104_ages.txt")
str(ages)
colnames(ages) <- ages[1,]
ages <- ages[-1,]
ages <- data.frame(apply(ages, 2, as.numeric)) #transform to numeric


#this is to transform to tidy format, calculate % and subset more common species
new <- forams %>% 
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
  mutate(total_sample = sum(count)) %>% 
  left_join(ages[c("depth", "min", "max", "mean")], by=c("depth")) %>%
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
  mutate(n_valves = ifelse(total_sample < 100, "low", "high")) #here label samples if nvalves >100 or higher
#drop_na()

#make it wide
core_counts_wide_forams <- core_counts_common %>%
  select(depth, mean, min,max, taxa, relative_abundance_percent) %>%
  rename(age_calyr=mean) %>%
  rename(upper_age=min) %>%
  rename(lower_age=max) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  arrange(depth) #sort by increasing time


is.na(core_counts_wide_forams$age_calyr)

# here assign manually age to depths that join did not work for some reason
core_counts_wide_forams[10,2] <- c(1581)
core_counts_wide_forams[10,3] <- c(1351)
core_counts_wide_forams[10,4] <- c(1785)
core_counts_wide_forams[75,2] <- c(8315)
core_counts_wide_forams[75,3] <- c(8089)
core_counts_wide_forams[75,4] <- c(8601)
core_counts_wide_forams[78,2] <- c(8601)
core_counts_wide_forams[78,3] <- c(8295)
core_counts_wide_forams[78,4] <- c(8970)
core_counts_wide_forams[80,2] <- c(8790)
core_counts_wide_forams[80,3] <- c(8445)
core_counts_wide_forams[80,4] <- c(9190)
core_counts_wide_forams[82,2] <- c(9026)
core_counts_wide_forams[82,3] <- c(8664)
core_counts_wide_forams[82,4] <- c(9431)
core_counts_wide_forams[83,2] <- c(9073)
core_counts_wide_forams[83,3] <- c(8700)
core_counts_wide_forams[83,4] <- c(9500)
core_counts_wide_forams[88,2] <- c(9738)
core_counts_wide_forams[88,3] <- c(9364)
core_counts_wide_forams[88,4] <- c(10202)
core_counts_wide_forams[91,2] <- c(9939)
core_counts_wide_forams[91,3] <- c(9586)
core_counts_wide_forams[91,4] <- c(10354)
core_counts_wide_forams[93,2] <- c(10073)
core_counts_wide_forams[93,3] <- c(9737)
core_counts_wide_forams[93,4] <- c(10483)
core_counts_wide_forams[95,2] <- c(10236)
core_counts_wide_forams[95,3] <- c(9935)
core_counts_wide_forams[95,4] <- c(10635)


plot.ts(diff(core_counts_wide_forams$upper_age))

## extract agedepth variables
agedepth <- core_counts_wide_forams[,names(core_counts_wide_forams) %in% c("depth", "age_calyr", "upper_age", "lower_age")]
forams <- core_counts_wide_forams[,!(names(core_counts_wide_forams) %in% c("upper_age", "lower_age", "depth", "age_calyr"))]

# Transform data to Hellinger form
forams[is.na(forams)] <- 0 #Replace NA (if any) by 0
forams_hell <- decostand(forams, method="hellinger")

# Run Principal Curves
forams_prc <- prcurve(forams, method = "ca", trace = TRUE, vary = TRUE, penalty = 1.4)

## Extract position on the curve
scrs_prc <- scores(forams_prc, display = "curve")

# Combine dataframe with ages and depths
foramsPrC <- cbind(agedepth, scrs_prc)

# Plot Pcurves with depth and ages
forams_plt_prc <- ggplot(foramsPrC, aes(x = age_calyr, y = PrC)) +
  geom_line() + geom_point() +
  labs(y = "PrC", x = "Age (cal yr BP)", title = "") +
  ggtitle("Forams PrC") +
  theme_bw()
forams_plt_prc


#Transform dataframe to include elapsedtime ("years mud slice") and standarized PrC scores
core <- transform(foramsPrC, Age = upper_age, negAge = - upper_age, elapsedTime = abs(upper_age - lower_age), 
                  rootPrC = sqrt(PrC), logPrC = log10(PrC+0.25))

## fit Gaussian location-scale model
mod <- gam(list(logPrC ~ s(negAge, k=30, bs="fs"), 
                      ~ elapsedTime + s(negAge)), #elapsedTime is the linear predictor
                 data = core, method = "REML", 
                 family = gaulss(link = list("identity", "logb")))


layout(matrix(1:4, ncol = 2))
gam.check(mod)

summary(mod)

## Predictions over range of both data sets
Nnew <- 500
elapsed <- 20 #evenly-spaced values over `Year`

createNewData <- function(Age, n, elapsed) {
  out <- data.frame(Age = seq(min(Age), max(Age), length.out = n))
  out <- transform(out, negAge = - Age, elapsedTime = elapsed)
  out
}
core <- core %>% drop_na()
newData <- with(core, createNewData(Age, Nnew, elapsed))

predGLSS <- function(model, newdata) {
  N <- nrow(newdata)
  p <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  fam <- family(model)
  mulink <- fam[["linfo"]][[1L]][["linkinv"]]
  sigmalink <- fam[["linfo"]][[2L]][["linkinv"]]
  res <- data.frame(Age = rep(newdata$Age, 2),
                    negAge = rep(newdata$negAge, 2),
                    term = rep(c("mu","sigma"), each = N),
                    fitted = c(mulink(p$fit[,1]), 1 / sigmalink(p$fit[,2])),
                    upper  = c(mulink(p$fit[,1] + (2 * p$se.fit[,1])),
                               1 / sigmalink(p$fit[,2] + (2 * p$se.fit[,2]))),
                    lower  = c(mulink(p$fit[,1] - (2 * p$se.fit[,1])),
                               1 / sigmalink(p$fit[,2] - (2 * p$se.fit[,2]))))
  res
}

predData <- predGLSS(mod, newData)

## plot
predForamsPlt <- ggplot(predData, aes(x = Age, y = fitted, group = term)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.1) +
  #geom_point(data= core, aes(x = Age, y = PrC), size=0.06) +
  geom_line() +
  facet_wrap(~ term, nrow = 2, labeller = label_parsed, scales = "free_y") +
  scale_x_reverse() +
  labs(y = "Fitted GAM-PrC", x = "Age (cal years BP)", title = "") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        strip.text.x = element_text(size=14))
predForamsPlt

ggsave("outputs/forams-GAM-PrC-variance.png", predForamsPlt, height = 8, width = 10)




