
rm(list=ls(all=TRUE))
dev.off()
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

##Load libraries for functions used
library(rbacon)
library(tidyverse)

## Build mixed curve for marine and brackish material
mix.ccurves(proportion = 0.5, cc1 = "3Col_intcal13.14C",
           cc2 = "3Col_marine20.14C", name = "mixed.14C", dirname = "Bacon_runs")

#Carlet core
Bacon("carlet", thick = 5, prob = 0.95, sep = ";",
      d.min = 74, d.max = 1620, d.by = 10, unit = "cm", cc=2, acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)

  ## mixed calibration curve
  Bacon(core = "carlet", thick = 5, prob = 0.95,
        d.min = 74, d.max = 1620, d.by = 10, unit = "cm", cc=4, cc4="mixed.14C", acc.shape = 2, 
        acc.mean = 30, mem.strength = 4, mem.mean = 0.7)

#S2 core
Bacon("S2_core_v2", thick = 20, prob = 0.95, sep = ";",
      d.by = 1, unit = "cm", acc.shape = 1.5, 
      acc.mean = 20, mem.strength = 10, mem.mean = 0.5, rotate.axes=TRUE)

add.dates(9550, 50, 3810)




# set parameters
thick<-20
d.by<-1
core.name="S2_core_v2"
acc.mean=20
acc.shape=1.5
mem.strength=10
mem.mean=0.5
# mem.strength * (1 - mem.mean) #should be smaller than 1

x <- Bacon(
  core.name,
  acc.mean=acc.mean,
  acc.shape=acc.shape,
  mem.mean=mem.mean,
  mem.strength=mem.strength,
  thick=thick,
  d.by=d.by,
  ask=FALSE,
  suggest=FALSE,
  plot.pdf=TRUE,
  depths.file=FALSE,
  normal=FALSE,
  rotate.axes=TRUE)


#functions to compute weighted mean and standard deviation (from SDMTools)
wt.mean <- function(x, wt) {
  s = which(is.finite(x*wt)); wt = wt[s]; x = x[s] #remove NA info
  sum(wt * x)/sum(wt)
}

wt.sd <- function(x, wt) {
  s = which(is.finite(x + wt)); wt = wt[s]; x = x[s] #remove NA info
  xbar = wt.mean(x,wt) #get the weighted mean
  sqrt(sum(wt *(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2))))
}


#extracting calibration object from info
calibration <- info$calib$probs

#every object in calibration is a matrix. The first column is the age, and the second column is the probability
calibration[[1]]

#getting ages and probs from the matrix
ages <- calibration[[1]][, 1]
probs <- calibration[[1]][, 2]

#calibration curve
plot(ages, probs, type = "l")

#weighted mean of first date
wt.mean(x = ages, wt = probs)

#weighted standard deviation of first date
wt.sd(x = ages, wt = probs)

#the important part

#weighted mean for every date
calibration.mean <- lapply(
  calibration,
  FUN = function(x) wt.mean(x = x[, 1], wt = x[, 2])
) %>%
  unlist()

#weighted sd for every date
calibration.sd <- lapply(
  calibration,
  FUN = function(x) wt.sd(x = x[, 1], wt = x[, 2])
) %>%
  unlist()

#to the dates data frame
df <- read.csv("Bacon_runs/S2_core_v2/S2_core_v2.csv",
  header = TRUE,
  sep = ";"
) %>%
  dplyr::mutate(
    calibration.mean = calibration.mean,
    calibration.sd = calibration.sd
  )

write.xlsx(df,"outputs/S2_calibrated.xlsx", row.names = FALSE)








## Examples from other lake sediment cores
#yahuarcocha core
Bacon(core="yahuarcocha", prob = 0.95, thick = 5,
      d.min = 1, d.max = 57, d.by = 1, unit = "cm", acc.shape = 1.5, acc.mean = 20, mem.strength = 4,
      mem.mean = 0.7, ssize = 5000, delta.R = 0)

#Pinan
Bacon(core="pinan", thick = 1, prob = 0.95,
      d.min = 1, d.max = 74, d.by = 1, unit = "cm", acc.shape = 1.5, acc.mean = 30, mem.strength = 4,
      mem.mean = 0.7)
  
  #Melina's code
  Bacon(core = "pinan", thick = 1, d.min = 1, d.max = 74,  acc.shape = 1.5, 
        acc.mean = c(50,100), slump = c(39,43,52,55))

  #new model 
  Bacon(core = "fondococha_new", thick = 1, d.min = 1, d.max = 75,  acc.shape = 1.5, 
        acc.mean = c(50,100), slump = c(39,40))
  
  
# Bacon Age Modelling to generate tables for MCEOF (code from Jon Tyler)
#get sample depths
#Llaviucu
proxy <- read.csv("/Volumes/xbenitogranell-data/0_project/data/cores/llaviucu_data.csv") %>%
  filter(depth<141)
depth <- proxy[,1]
original_age <- read.table("Bacon_runs/llaviucu/llaviucu_28_ages.txt", header=TRUE)[,5]

#Yahuarcocha
proxy <- read.csv("/Volumes/xbenitogranell-data/0_project/data/cores/yahuarcocha_up.csv")
depth <- proxy[,3]
original_age <- read.table("Bacon_runs/yahuarcocha/yahuarcocha_12_ages.txt", header = TRUE)

#Pinan
proxy <- read.csv("/Volumes/xbenitogranell-data/0_project/data/cores/pinan_up.csv")
depth <- proxy[,2]

#Fondococha
proxy <- read.csv("/Volumes/xbenitogranell-data/0_project/data/cores/fondococha_up.csv") #counts
depth <- proxy[,1]


#original_age <- proxy[,5]

# define how many iterations to do 
nITS <- 10000

# create a matrix to store results
ages <- matrix(NA, length(depth), nITS)

# get nITS age-depth iterations and save them as a .csv file
for(i in 1:length(depth)){
  DEPTH <- depth[i]
  ages[i, ] <- Bacon.Age.d(DEPTH)[1:nITS]										
}
write.csv(ages, "llaviucu_ages.csv")
write.csv(ages, "yahuarcocha_ages.csv")
write.csv(ages, "pinan_ages.csv")
write.csv(ages, "fondococha_ages.csv")

write.csv(diat, "llaviucu_proxy.csv")
write.csv(diat, "yahuarcocha_proxy.csv")
write.csv(diat, "pinan_proxy.csv")
write.csv(diat, "fondococha_proxy.csv")

lines(depth, original_age, col="orange", lwd=2)

#############



#Age-depth modeling Clam-Stars (Flantua)
setwd("/nfs/xbenitogranell-data/0_project/R codes/Clam_stars")
setwd("/Volumes/xbenitogranell-data/0_project/R codes/Clam_stars")

source("Clam.R")
source("StarClassification_AgeModels.R")

#pinan
Agemodel.stars("pinan", cc=1) 

#yahuarcocha
Agemodel.stars()
