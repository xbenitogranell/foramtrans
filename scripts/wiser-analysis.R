## Analyse Hel's WISER ordination data

## load packages
library("analogue")
library("mgcv")

## Load the data
source("core-data-processing.R")

## Load functions
source("functions.R")

## need time tracks for each site
doTimeTrack <- function(core, train, scale = TRUE,
                        scaling = 3) {
    DROP <- which(colnames(core) %in% c("sites","depths","dates"))
    timetrack(train, core[, -DROP], method = "rda", scale = scale,
              scaling = scaling)
}

## split up by sites
coresSpl <- split(cores, cores$sites)
ttracks <- lapply(coresSpl, doTimeTrack, train = deep,
                  scale = FALSE, scaling = 1)


## fit the surface once - it is the same for each deep or each shallow lake
surfDeep <- ordisurf(ttracks[[1]]$ord, deepTP[, "Logtp"],
                     method = "REML", scaling = ttracks[[1]]$scaling,
                     select = TRUE, plot = FALSE, knots = 20)

## plot timtracks
SURF <- NULL # surfDeep                  # or NULL to turn off
tt.lwd <- 1.5                     # lwd for plots, smaller for PDF/EPS
for(i in seq_along(ttracks)) {
    fname <- paste0("./fig/eps/deep-", tolower(names(ttracks)[i]),
                    "-timetrack.eps")
    postscript(fname, pointsize = 12, height = 8, width = 8, onefile = FALSE,
               horizontal = FALSE)
    op <- par(mar = c(5,4,1,1) + 0.1)
    ttplot(ttracks[[i]], surf = SURF, lwd = tt.lwd)   # ttplot load from functions.R
    par(op)
    dev.off()
}

## Shallow lakes ------------------------------------------------------
## split up by sites
shcoresSpl <- split(shcores, shcores$sites)
shttracks <- lapply(shcoresSpl, doTimeTrack, train = shallow,
                    scale = FALSE, scaling = 1)

## fit the surface once
surfShallow <- ordisurf(shttracks[[1]]$ord, shallowTP[, "Log10tp"],
                        method = "REML", scaling = shttracks[[1]]$scaling,
                        select = TRUE, plot = FALSE, knots = 20)

## plot timtracks
for(i in seq_along(shttracks)) {
    fname <- paste0("./fig/eps/shallow-", tolower(names(shttracks)[i]),
                    "-timetrack.eps")
    postscript(fname, pointsize = 12, height = 8, width = 8, onefile = FALSE,
               horizontal = FALSE)
    op <- par(mar = c(5,4,1,1) + 0.1)
    ttplot(shttracks[[i]], surf = SURF, lwd = tt.lwd)   # ttplot load from functions.R
    par(op)
    dev.off()
}

## plot the species codes
plot(shttracks[[1]]$ord, display = "species", scaling = 3,
     xlim = c(-1,0.8), ylim = c(-0.6,0.8))

## which taxa are abundant or common?
postscript("./fig/eps/shallow-species-plot.eps", height = 5, width = 10,
           onefile = FALSE, horizontal = FALSE, pointsize = 10, paper = "special")
layout(matrix(1:2, ncol = 2))
op <- par(mar = c(5,4,1,1) + 0.1)
set.seed(24)
scrs <- scores(shttracks[[1]]$ord, display = "species", scaling = 3)
take <- colnames(shallow) %in% rownames(scrs)
TAXA <- which(colSums(shallow[,take] > 0) > 5 & (apply(shallow[,take]^2, 2, max) > 5) &
                  ((scrs[,1] <= -0.5 | scrs[,1] >= 0.5) | (scrs[,2] <= -0.5 | scrs[,2] >= 0.5))
              )
plot(shttracks[[1]]$ord, display = "species", scaling = 3, type = "n")
rect(-1, -0.6, 0.8, 0.8, border = "grey", col = NA)
points(shttracks[[1]]$ord, display = "species", cex = 0.5, scaling = 3)
ordipointlabel(shttracks[[1]]$ord, display = "species", scaling = 3,
               select = TAXA, cex = 0.5, add = TRUE)
TAXA <- which(colSums(shallow[,take] > 0) > 5 & (apply(shallow[,take]^2, 2, max) > 5))
plot(shttracks[[1]]$ord, display = "species", scaling = 3, type = "n",
     xlim = c(-1,0.8), ylim = c(-0.6,0.8))
rect(-1, -0.6, 0.8, 0.8, border = "grey", col = NA)
points(shttracks[[1]]$ord, display = "species", cex = 0.5, scaling = 3)
ordipointlabel(shttracks[[1]]$ord, display = "species", scaling = 3,
               xlim = c(-1,0.8), ylim = c(-0.6,0.8), select = TAXA,
               cex = 0.5, add = TRUE)
par(op)
layout(1)
dev.off()

plot(ttracks[[1]]$ord, display = "species", scaling = ttracks[[1]]$scaling)

postscript("./fig/eps/deep-species-plot.eps", height = 5, width = 10,
           onefile = FALSE, horizontal = FALSE, pointsize = 10, paper = "special")
layout(matrix(1:2, ncol = 2))
op <- par(mar = c(5,4,1,1) + 0.1)
set.seed(24)
scrs <- scores(ttracks[[1]]$ord, display = "species", scaling = 3)
take <- colnames(deep) %in% rownames(scrs)
TAXA <- which(colSums(deep[,take] > 0) > 2 & (apply(deep[,take]^2, 2, max) > 2) &
                  ((scrs[,1] <= -0.35 | scrs[,1] >= 0.5) | (scrs[,2] <= -0.35 | scrs[,2] >= 0.25)))
plot(ttracks[[1]]$ord, display = "species", scaling = 3, type = "n")
rect(-0.35, -0.35, 0.5, 0.25, border = "grey", col = NA)
points(ttracks[[1]]$ord, display = "species", cex = 0.5, scaling = 3)
ordipointlabel(ttracks[[1]]$ord, display = "species", scaling = 3,
               select = TAXA, cex = 0.5, add = TRUE)
TAXA <- which(colSums(deep[,take] > 0) > 2 & (apply(deep[,take]^2, 2, max) > 2))
plot(ttracks[[1]]$ord, display = "species", scaling = 3, type = "n",
     xlim = c(-0.35, 0.5), ylim = c(-0.35, 0.25))
rect(-0.35, -0.35, 0.5, 0.25, border = "grey", col = NA)
points(ttracks[[1]]$ord, display = "species", cex = 0.5, scaling = 3)
ordipointlabel(ttracks[[1]]$ord, display = "species", scaling = 3,
               xlim = c(-1,1.75), ylim = c(-0.75,0.75),
               select = TAXA, cex = 0.5, add = TRUE)
par(op)
layout(1)
dev.off()
