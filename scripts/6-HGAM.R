
#Clear workspace
rm(list=ls(all=TRUE))

##loading libraries for functions used
library(analogue) 
library(rioja) 
library(plyr) #allow to join dataframes by common column
library(dplyr) #allow to summarise variables and manipulate multiple dataframes
library(ggplot2) #to make nice plots
library(tidyverse)
library(tidypaleo) #(https://fishandwhistle.net/post/2018/stratigraphic-diagrams-with-tidypaleo-ggplot2/)
library(cluster)
library(mgcv)

## Read S2 foraminifera records
forams <- read.csv("datasets/S2/S2_counts.csv", sep=";")[-1] %>%
  mutate(depth=as.numeric(gsub(",", ".", gsub("\\.", "", depth))),
         depth=depth*100) #replace commas with dots for decimals
str(forams)

# Read in forams taxa names groups
#changes <- read.csv("datasets/S2/nms_taxa_groups.csv")

# Read nms list to standardize taxa names
nms <- read.csv("outputs/nms_old_new.csv") %>%
  select(c(1:3))

# Read in age-depth model
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

# 
# core_counts_wide_forams[c(18,19),2] <- 1581
# core_counts_wide_forams[c(203,204,205),2] <- 8315
# core_counts_wide_forams[c(209),2] <- 8601
# core_counts_wide_forams[c(212),2] <- 8790
# core_counts_wide_forams[c(214),2] <- 9026
# core_counts_wide_forams[c(215,216),2] <- 9073
# core_counts_wide_forams[c(228,229),2] <- 9738
# core_counts_wide_forams[c(235,236),2] <- 9939
# core_counts_wide_forams[c(238,239),2] <- 10073
# core_counts_wide_forams[c(241,242),2] <- 10236
# 
# 
# 
# core_counts_wide_forams[c(203,204,205),2] <- 8315
# core_counts_wide_forams[c(209),2] <- 8601
# core_counts_wide_forams[c(212),2] <- 8790
# core_counts_wide_forams[c(214),2] <- 9026
# core_counts_wide_forams[c(215,216),2] <- 9073
# core_counts_wide_forams[c(228,229),2] <- 9738
# core_counts_wide_forams[c(235,236),2] <- 9939
# core_counts_wide_forams[c(238,239),2] <- 10073
# core_counts_wide_forams[c(241,242),2] <- 10236
# 
# core_counts_wide_forams[c(18,19),c(3,4)] <- c(1351,1785)
# core_counts_wide_forams[c(203),c(3)] <- c(8089)
# core_counts_wide_forams[c(203),c(4)] <- c(8601)
# core_counts_wide_forams[c(204),c(3)] <- c(8089)
# core_counts_wide_forams[c(204),c(4)] <- c(8601)
# core_counts_wide_forams[c(205),c(3)] <- c(8089)
# core_counts_wide_forams[c(205),c(4)] <- c(8601)
# core_counts_wide_forams[c(209),c(3)] <- c(8295)
# core_counts_wide_forams[c(209),c(4)] <- c(8970)
# core_counts_wide_forams[c(212),c(3)] <- c(8445)
# core_counts_wide_forams[c(212),c(4)] <- c(9190)
# core_counts_wide_forams[c(214),c(3)] <- c(8664)
# core_counts_wide_forams[c(214),c(4)] <- c(9431)
# core_counts_wide_forams[c(215,216),c(3,4)] <- c(8700,9500)
# core_counts_wide_forams[c(228,229),c(3,4)] <- c(9364,10202)
# core_counts_wide_forams[c(235,236),c(3,4)] <- c(9586,10354)
# core_counts_wide_forams[c(238,239),c(3,4)] <- c(9737,10483)
# core_counts_wide_forams[c(241,242),c(3,4)] <- c(9935,10635)

# Read Carlet and St Jaume cores
carlet_counts_wide <- read.csv("datasets/Carlet/carlet_counts_wide.csv")[-1]
stjaume_counts_wide <- read.csv("datasets/Sant Jaume/stjaume_counts_wide.csv")[-1]

core_counts_wide_forams <- stjaume_counts_wide

#Prepare the data
agedepth <- core_counts_wide_forams[, names(core_counts_wide_forams) %in% c("depth", "age_calyr", "upper_age","lower_age")] 
forams <- core_counts_wide_forams[, !names(core_counts_wide_forams) %in% c("depth", "age_calyr", "upper_age","lower_age")]
forams[is.na(forams)] <- 0

#Select most common species 
criteria <- 0.2 #% of the total samples

n.occur <- apply(forams>0, 2, sum)
forams_red <- forams[, n.occur > (dim(forams)[1])*criteria] #
forams <- cbind(agedepth, forams_red)

#this is to transform to tidy format, calculate % and subset more common species
forams_data <- forams %>% 
  gather(key = taxa, value = count, -depth, -age_calyr, -upper_age, -lower_age) %>%
  group_by(depth) %>%
  mutate(total_sample = sum(count)) %>% 
  filter(!total_sample == "0") %>% #this is to remove empty samples
  mutate(log_total_counts = log10(total_sample+1)) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  mutate(negAge = -age_calyr) %>%
  mutate(elapsedTime = round(abs(upper_age - lower_age),0)) %>%
  ungroup() %>%
  filter(!is.na(elapsedTime)) %>%
  mutate(spp = factor(taxa)) 

levels(forams_data$spp)

#model S HGAM : similar smootheness between groups (spp) without global smooth 
set.seed(10) #set a seed so this is repeatable
gam_S <- gam(count ~ s(negAge, spp, k=30, bs="fs") + offset(log_total_counts),
                    weights = elapsedTime/mean(elapsedTime),
                    data=forams_data, family = nb, 
                    method = "REML")

summary(gam_S)
gam.check(gam_S)
#draw(gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
gam_I<- gam(count ~ s(negAge, by=spp, k=30, bs="fs") +
                     s(spp, bs="re") + offset(log_total_counts),
                   weights = elapsedTime/mean(elapsedTime),
                   data=forams_data, family = nb, #places notes at the deciles of sample ages
                   method = "REML")

gam.check(gam_I)
draw(gam_I)

#Compare different model fits using AIC
AIC_table <- AIC(gam_S, gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

#Create synthetic data to predict over a range of ages
forams_plot_data <- with(forams_data, as_tibble(expand.grid(negAge = seq(min(forams_data$negAge), max(forams_data$negAge)),
                                                        spp = factor(levels(forams_data$spp)),
                                                        log_total_counts = mean(log_total_counts))))

modS_fit <- predict(gam_S, newdata = forams_plot_data, se.fit = TRUE)
modI_fit <- predict(gam_I, newdata = forams_plot_data, se.fit = TRUE)

#non-shared trends
forams_plot_data$modS_fit <- as.numeric(modS_fit$fit)
forams_plot_data$modI_fit <- as.numeric(modI_fit$fit)

# comparing non-shared trends
forams_plot_data <- gather(forams_plot_data, key=model, value=fit, modS_fit, modI_fit)

forams_plot_data <- mutate(forams_plot_data, se= c(as.numeric(modS_fit$se.fit),
                                               as.numeric(modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

#Plot the model output for non-shared trends, with means plus standard deviations for each model.
forams_plot_model_labels <- paste("Model", c("S", "I"))
forams_plot_model_labels <- factor(forams_plot_model_labels, levels = forams_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

forams_plot <- ggplot(forams_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= forams_data, aes(x = negAge, y = count), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Absolute counts", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = forams_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = forams_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

forams_plot

## sAVE THE PLOT
ggsave("outptus/forams-HGAM.png", forams_plot, height = 8, width = 10)

## save model results for later use
write.csv(forams_plot_data, "outputs/stjaume-forams-HGAMs-fitted-values.csv", row.names = FALSE)


##Derivatives and posterior distribution simulation

# Eric's code
#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 180

years <- seq(min(forams_plot_data$negAge),
             max(forams_plot_data$negAge),
             length.out = n_length)

diff(years)

#model <- gam_S
model <- gam_I

#pred <- modS_fit
pred <- modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(forams_plot_data$spp),
                         negAge = seq(min(forams_plot_data$negAge),
                                      max(forams_plot_data$negAge),
                                      length.out = n_length),
                         log_total_counts=0)

map_pred_sims <- predict(model,
                         confint_sims,
                         type = "lpmatrix") %*% random_coefs %>%
  as_data_frame() %>%
  bind_cols(confint_sims)%>%
  gather(key = simulation, value = pred, -negAge, -log_total_counts,-spp)


#specifying the step size for numerical derivative calculations
delta = 0.01

#calculating the predicted value for the current year plus delta
step_ahead_fits = confint_sims %>%
  mutate(negAge = negAge+delta)%>%
  predict(model, 
          ., type = "lpmatrix") %*% random_coefs 


#calculating the predicted value for the current year minus delta
step_behind_fits = confint_sims %>%
  mutate(negAge = negAge-delta)%>%
  predict(model,
          ., type = "lpmatrix") %*% random_coefs 


#using the predicted values for year plus and minus delta to calculate
#derivatives for each species for each simulation
derivs <- calc_1st_deriv(step_behind_fits,step_ahead_fits,delta = delta)%>%
  as_data_frame()%>%
  bind_cols(confint_sims)%>%
  gather(key = simulation,value = deriv, -spp,-negAge, -log_total_counts)

#Creating summaries of derivatives for each simulation for each year
deriv_summaries <- derivs %>%
  group_by(negAge,simulation)%>%
  summarize(deriv_mean = mean(deriv),
            deriv_sd = sd(deriv))%>%
  group_by(negAge)%>% #turning derivative summaries into 95% confidence intervals
  select(-simulation)%>%
  summarize_all(.funs = list(lower = ~quantile(.,probs = 0.025),
                             upper = ~quantile(.,probs = 0.975),
                             med   = ~quantile(.,probs = 0.5)))

#Plotting mean rate of change plus the 95% CI
mean_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_mean_med, 
             ymin = deriv_mean_lower,
             ymax = deriv_mean_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
mean_plot 

#Plotting standard deviation of rate of change plus the 95% CI
sd_plot <- deriv_summaries %>%
  ggplot(aes(x = negAge, 
             y = deriv_sd_med, 
             ymin=deriv_sd_lower,
             ymax=deriv_sd_upper))+
  geom_ribbon(fill="grey")+
  geom_line()+
  geom_hline(yintercept = 0, linetype=2) +
  scale_y_continuous("")+
  xlab("Cal years BP") +
  theme_bw()
sd_plot


## save derivative summaries for later use
write.csv(deriv_summaries, "outputs/stjaume-forams-derivatives.csv", row.names = FALSE)
