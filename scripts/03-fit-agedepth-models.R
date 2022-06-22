

rm(list=ls(all=TRUE))
dev.off()
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

##Load libraries for functions used
library(rbacon)
library(tidyverse)

## Build mixed curve for marine and brackish material
mix.curves(proportion = 0.5, cc1 = "3Col_intcal13.14C",
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
Bacon("S2_core", thick = 20, prob = 0.95, sep = ";",
      d.min = 1800, d.max = 3810, d.by = 10, unit = "cm", cc=2, acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)

add.dates(9550, 50, 3810)



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
