

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


Bacon("carlet", thick = 5, prob = 0.95, sep = ";",
      d.min = 74, d.max = 1620, d.by = 10, unit = "cm", cc=2, acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)

## mixed calibration curve
Bacon(core = "carlet", thick = 5, prob = 0.95,
      d.min = 74, d.max = 1620, d.by = 10, unit = "cm", cc=4, cc4="mixed.14C", acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)



## Examples from lake sediment cores
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
  
  

#Titicaca (Weide et al 2017) CHUA3 core
Bacon("titicaca", thick = 4, prob = 0.95,
      d.min = 0, d.max = 169, d.by = 1, unit = "cm", cc=4, cc4="mixed.14C", hiatus.depths = c(47, 54, 73), acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)

#Titicaca (Weide et al 2017) LT01-3A core
Bacon("titicaca-LT01-3A", thick = 1, prob = 0.95, cc=3,
      unit = "cm", d.min = 117, d.max=4929, hiatus.depths = c(290, 344, 390),
      acc.shape = 1.5, acc.mean = 15, mem.strength = 4, mem.mean = 0.7)


#Umayo (Ekhdal et al 2008)
Bacon("umayo", thick = 2, prob = 0.95,
      d.min = 1, d.max = 678, d.by = 0.1, unit = "cm", acc.shape = 2, 
      acc.mean = 30, mem.strength = 4, mem.mean = 0.7)


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
