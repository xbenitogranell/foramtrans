# Function to transform relative abundance to absolute counts

## Read in forams counts
forams <- read.csv("datasets/Sant Jaume/stjaume_forams_RA.csv", sep=",")
forams <- read.csv("datasets/Carlet/carlet_forams_RA.csv", sep=",")

mydata_num <- data.frame(apply(forams[,3:ncol(forams)], 2, as.numeric)) 
mydata_num[is.na(mydata_num)] <- 0
str(mydata_num)
colSums(mydata_num)
rowSums(mydata_num)

prop <- apply(mydata_num, 1, FUN = function(x) {min(x[x > 0])})
prop
fc <- 1/prop
fc

counts <- data.frame(apply(mydata_num, 2, FUN = function(x) {round(x*fc,0)})) #when >1 observations
rowSums(counts)
str(counts)

stjaume_forams <- data.frame(depth=forams[,2],counts)
carlet_forams <- data.frame(depth=forams[,2], counts)
write.csv(stjaume_forams, "datasets/Sant Jaume/stjaume_forams_counts.csv")
write.csv(carlet_forams, "datasets/Carlet/carlet_forams_counts.csv")
