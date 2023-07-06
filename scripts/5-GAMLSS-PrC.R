

#function to perform GAM Gausian location and scale
fitGamLS <- function(i, dat, kk) {
  cdata <- dat[[i]]
  k <- kk[i]
  fit <- gam(list(P_B ~ s(negAge, k = k, bs="ad"), ~ elapsedTime + s(negAge)),
             data = cdata, method = "REML", 
             family = gaulss(link = list("identity", "logb")))
  
}

k <- c(rep(20,8))
fits <- lapply(seq_along(coreDataSplit), fitGamLS, dat = coreDataSplit, kk = k)
names(fits) <- names(coreDataSplit)


# Function to plot GAM smooths
plotSmooths <- function(i, fits, ...) {
  plot(fits[[i]], ..., main = names(fits)[i])
}

layout(matrix(1:16, ncol = 4))
op <- par(mar = c(3,4,1,0.4) + 0.1)
lapply(seq_along(fits), plotSmooths, fits = fits, xlab = "", title=names(fits), residuals = TRUE,
       pch = 1, cex = 0.4)
par(op)
layout(1)


# Make predictions
Nnew <- 500
elapsed <- 20 #evenly-spaced values over `Year`

createNewData <- function(Age, n, elapsed) {
  out <- data.frame(Age = seq(min(Age), max(Age), length.out = n))
  out <- transform(out, negAge = - Age, elapsedTime = elapsed)
  out
}

core <- coresRatios
newDiat <- with(core, createNewData(Age, Nnew, elapsed))


predData <- lapply(seq_along(coreDataSplit), createNewData(Age, Nnew, elapsed))


makePredData <- function(data, n = 200) {
  data.frame(negAge = seq(min(data$negAge), max(data$negAge), length.out = n))
}

predData <- lapply(coreDataSplit, makePredData)


predGLSS <- function(i, models, newdata) {
  N <- nrow(newdata)
  p <- predict(models[[i]], newdata = newdata[[i]], type = "link", se.fit = TRUE)
  fam <- family(models[[i]])
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

predDiat <- predGLSS(fits, newDiat)

predTrends <- lapply(seq_along(fits), predGLSS,
                     models = fits, newdata = newDiat)


modelPreds <- function(i, models, newdata, se.fit = TRUE) {
  predict(models[[i]], newdata = newdata[[i]], se.fit = se.fit)
}

predTrends <- lapply(seq_along(fits), modelPreds,
                     models = fits, newdata = predData)




## GAMLS   
## plot the data
core <- coresPrC %>% 
  mutate(elapsedTime = abs(upper_age - lower_age)) %>%
  mutate(rootPrC = sqrt(PrC)) %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  filter(AgeCE >= 0) %>%
  filter(!upper_age == 0)

coreDataSplit <- split(core, core$Lake)


#Load data

##Subsetting
select <- c("llaviucu")
select <- c("yahuarcocha")
select <- c("titicaca")
select <- c("umayo")
select <- c("pinan")
select <- c("triumfo")

core <- coresPrC %>% 
  filter(str_detect(Lake, select)) %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  filter(AgeCE >= "0") 



## plot the data
diatomPlt <- ggplot(core, aes(x = Age, y = PrC)) +
  geom_point() +
  geom_line() +
  scale_x_reverse()
diatomPlt


#Transform dataframe to include elapsedtime ("years mud slice") and standarized PrC scores
core <- transform(core, Age = upper_age, negAge = - upper_age, elapsedTime = abs(upper_age - lower_age), 
                  rootPrC = sqrt(PrC), logPrC = log10(PrC+0.25))

## fit Gaussian location-scale model
modDiatom <- gam(list(PrC ~ s(negAge, k=30, bs="fs"), 
                      ~ elapsedTime + s(negAge)), #elapsedTime is the linear predictor
                 data = core, method = "REML", 
                 family = gaulss(link = list("identity", "logb")))


layout(matrix(1:16, ncol = 4))
gam.check(modDiatom)

summary(modDiatom)


## Predictions over range of both data sets
Nnew <- 500
elapsed <- 20 #evenly-spaced values over `Year`

createNewData <- function(Age, n, elapsed) {
  out <- data.frame(Age = seq(min(Age), max(Age), length.out = n))
  out <- transform(out, negAge = - Age, elapsedTime = elapsed)
  out
}
newDiat <- with(core, createNewData(Age, Nnew, elapsed))

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

predDiat <- predGLSS(modDiatom, newDiat)



## plot
predDiatPlt <- ggplot(predDiat, aes(x = Age, y = fitted, group = term)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.1) +
  #geom_point(data= core, aes(x = Age, y = PrC), size=0.06) +
  geom_line() +
  facet_wrap(~ term, nrow = 2, labeller = label_parsed, scales = "free_y") +
  scale_x_reverse() +
  labs(y = "Fitted", x = "Age (cal years BP)", title = "")
predDiatPlt

