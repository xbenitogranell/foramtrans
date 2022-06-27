## PCA-Nipals XRF
library(ade4)

# Read in XRF data
S2_geochem_data <- read.csv("datasets/S2/S2_XRF.csv") %>%
  dplyr::select(-LE) %>%
  dplyr::rename(depth=ï..depth )

str(S2_geochem_data)

# Prepare data
PCA.data <- data.frame(scale(S2_geochem_data[,2:ncol(S2_geochem_data)])) #scale variables

#Perform PCA with NIPALS algorithm
PCA.nipals <- nipals(PCA.data, nf = 2, rec = FALSE, niter = 1000, tol = 1e-09)

#Save summary matrix results  
PCA.summary <- NULL
PCA.summary <- matrix(data = NA, nrow = 6, ncol = 3, byrow = FALSE, dimnames = NULL)
colnames(PCA.summary) <- c("Parameter", "Value", "Explained variance")

PCA.summary[1, 1] <- c("Number of extracted factors")
PCA.summary[1, 2] <- PCA.nipals$nf

PCA.summary[2, 1] <- c("Number of variables")
PCA.summary[2, 2] <- nrow(PCA.nipals$co)

PCA.summary[3, 1] <- c("Eigenvalue 1")
PCA.summary[3, 2] <- PCA.nipals$eig[1]
PCA.summary[3, 3] <- (PCA.nipals$eig[1] / nrow(PCA.nipals$co)) * 100

PCA.summary[4, 1] <- c("Eigenvalue 2")
PCA.summary[4, 2] <- PCA.nipals$eig[2]
PCA.summary[4, 3] <- (PCA.nipals$eig[2] / nrow(PCA.nipals$co)) * 100

PCA.summary[5, 1] <- c("Number of iterations axis 1")
PCA.summary[5, 2] <- PCA.nipals$nb[1]

PCA.summary[6, 1] <- c("Number of iterations axis 2")
PCA.summary[6, 2] <- PCA.nipals$nb[2]

#Save the column coordinates (Component matrix) #ncol= n variables 
Component.matrix <- NULL
Component.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 2, byrow = FALSE, dimnames = NULL)
Component.matrix[,1] <- rownames(PCA.nipals$co)
Component.matrix[,2] <- (((PCA.nipals$co[,1] ^ 2) * PCA.nipals$eig[1]) / PCA.nipals$eig[1]) + (((PCA.nipals$co[,2] ^ 2) * PCA.nipals$eig[2]) / PCA.nipals$eig[2])
Component.matrix <- cbind(Component.matrix[,1], PCA.nipals$co, Component.matrix[,2])
colnames(Component.matrix) <- c("Variable", "Component 1", "Component 2", "Power Importance")

#Save the column normed scores (Component scores coefficient matrix)
Component.coefficient.matrix <- NULL
Component.coefficient.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 1, byrow = FALSE, dimnames = NULL)
Component.coefficient.matrix[,1] <- rownames(PCA.nipals$c1)
Component.coefficient.matrix <- cbind(Component.coefficient.matrix, PCA.nipals$c1)
colnames(Component.coefficient.matrix) <- c("Variable", "Component 1", "Component 2")

#Save the row coordinates (Factor Scores)
colnames(PCA.nipals$li) <- c("PCA1", "PCA2")
Factor.scores <- data.frame(cbind(PCA.data, PCA.nipals$li))
Factor.scores$depth <- S2_geochem_data$depth




#Create data frame with site scores and regions
#region <- site.time[1:21,] #lake regions of the modern diatom dataset
PCA.scores <- data.frame(component1=Factor.scores$PCA1, component2=Factor.scores$PCA2)

#extract factor scores (environmental variables)
comp1 <- as.numeric(Component.coefficient.matrix[,2])
comp2 <- as.numeric(Component.coefficient.matrix[,3])

#Plot PCA site labels (=lakes)
pca_tbl <- mutate(data.frame(PCA.scores))

pca_plt <- ggplot(pca_tbl, aes(component1,component2)) + 
  xlab("PCA1") + ylab("PCA2") +
  coord_fixed() +
  geom_point() +
  #geom_text_repel(colour="black", size=3) +
  geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  theme_classic()
pca_plt

#Plot PCA species (=environmental variables)  
variables_tbl <- mutate(data.frame(cbind(comp1, comp2)), varlbls = as.character(Component.coefficient.matrix[,1]))
pca_variables_plt <- ggplot(variables_tbl, aes(comp1,comp2, label=varlbls)) + 
  xlab("PCA1") + ylab("PCA2") +
  geom_point() +
  geom_text_repel(colour="black", size=3) +
  geom_segment(data=variables_tbl, aes(x = 0, y = 0, xend = comp1*0.9, yend = comp2*0.9), arrow = arrow(length = unit(1/2, 'picas')), color = "grey30") +
  geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  theme_classic()
pca_variables_plt

### Plot CCAs (ggplot)
## Run PCA with vegan
#check rows with NA and remove
row.has.na <- apply(PCA.data, 1, function(x){any(is.na(x))})
sum(row.has.na)
PCA.data <- na.omit(PCA.data)

# Run PCA with vegan
modPCA <- rda(PCA.data, scale=TRUE)
plot(modPCA, scale=3)

# Function to extract % of explained variability for CCA axes 1 and 2
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

(labs_PCA<- axis.expl(modPCA))

# Fortify the ordinations for ggploting
ford <- fortify(modPCA, axes = 1:2)  # fortify the ordination
take <- c('PC1', 'PC2')  # which columns contain the scores we want
sites <- subset(ford, Score == 'sites')  # take only biplot arrow scores
arrows <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'species'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

ford[locate, ]

#plot
PCA_XRF <- ggplot() +
  theme_bw()+
  geom_point(data = subset(ford, Score == 'sites'),
             mapping = aes(x = PC1, y = PC2)) +
  # geom_text(data = subset(ford, Score == 'sites'),
  #           mapping = aes(label=Label, x = CCA1, y = CCA2),size=3) + 
  geom_segment(data=arrows,
               mapping = aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.01, "npc")), colour="blue")+
  geom_text(data = arrows, # crudely push labels away arrow heads
            mapping = aes(label = Label, x = PC1 * 1.1, y = PC2 * 1.1), colour="blue") +
  xlab(paste0(names(labs_PCA[1]), " (", sprintf("%.1f", labs_PCA[1]), "%)"))+
  ylab(paste0(names(labs_PCA[2]), " (", sprintf("%.1f", labs_PCA[2]), "%)"))+
  coord_fixed()
PCA_XRF

ggsave("outputs/PCA_XRF.png", PCA_XRF, height = 6, width = 10)


## plot downcore PCA-NIpals factor scores
PCA_components_long <- gather(data=Factor.scores %>% select(c("PCA1", "PCA2","depth")), key = param, value = value, -depth) 

PCA_components_plt <- ggplot(PCA_components_long, aes(x = depth, y = value)) +
  geom_line() +
  geom_point() +
  scale_y_reverse() +
  scale_x_reverse() +
  facet_grid(~param) +
  coord_flip() +
  #geom_smooth() +
  labs(x = "core depth (m)", y = NULL) +
  ggtitle("")

