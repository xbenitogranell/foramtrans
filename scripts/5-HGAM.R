
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


## Calculate relative abundance
# sample_info <- forams[, names(forams) %in% c("sample_id", "depth")]
# forams <- forams[, !names(forams) %in% c("sample_id", "depth")]
# forams[is.na(forams)] <- 0
# 
# # Transform to relative abundance
# total <- apply(forams, 1, sum)
# forams_red <- forams[total>0, ] 
# 
# forams <- forams/total*100
# #
# # ##Remove rare species
# 
# abund <- apply(forams, 2, max)
# n.occur <- apply(forams>0, 2, sum)
# forams <- forams[, n.occur>1 & abund>2] #more than 2% of RA and present in >1 sample


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
  select(depth, mean, min,max, taxa, assemblage, count) %>%
  rename(age_calyr=mean) %>%
  rename(upper_age=min) %>%
  rename(lower_age=max) %>%
  spread(key = taxa, value = count) %>%
  arrange(depth) #sort by increasing time


is.na(core_counts_wide_forams$age_calyr)

# here assign manually age to depths that join did not work for some reason
core_counts_wide_forams[c(18,19),2] <- 1581
core_counts_wide_forams[c(203,204,205),2] <- 8315
core_counts_wide_forams[c(209),2] <- 8601
core_counts_wide_forams[c(212),2] <- 8790
core_counts_wide_forams[c(214),2] <- 9026
core_counts_wide_forams[c(215,216),2] <- 9073
core_counts_wide_forams[c(228,229),2] <- 9738
core_counts_wide_forams[c(235,236),2] <- 9939
core_counts_wide_forams[c(238,239),2] <- 10073
core_counts_wide_forams[c(241,242),2] <- 10236


#Prepare the data
agedepth <- core_counts_wide_forams[, names(core_counts_wide_forams) %in% c("depth", "age_calyr", "assemblage")] 
forams <- core_counts_wide_forams[, !names(core_counts_wide_forams) %in% c("depth", "age_calyr", "assemblage")]
forams[is.na(forams)] <- 0

#Select most common species 
criteria <- 0.1 #% of the total samples

n.occur <- apply(forams>0, 2, sum)
forams_red <- forams[, n.occur > (dim(forams)[1])*criteria] #
forams <- cbind(agedepth, forams_red)

#this is to transform to tidy format, calculate % and subset more common species
diat_data <- forams %>% 
  gather(key = taxa, value = count, -depth, -assemblage, -age_calyr) %>%
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

levels(diat_data$spp)

#model S HGAM : similar smootheness between groups (spp) without global smooth 
set.seed(10) #set a seed so this is repeatable
diatom_gam_S <- gam(count ~ s(negAge, spp, k=20, bs="fs") + offset(log_total_counts),
                    weights = elapsedTime/mean(elapsedTime),
                    data=diat_data, family = nb, #places notes at the deciles of sample ages
                    method = "REML")

summary(diatom_gam_S)
gam.check(diatom_gam_S)
draw(diatom_gam_S)

#model I HGAM: different smootheness for each taxa without global smooth
diatom_gam_I<- gam(count ~ s(negAge, by=spp, k=20, bs="fs") +
                     s(spp, bs="re") + offset(log_total_counts),
                   weights = elapsedTime/mean(elapsedTime),
                   data=diat_data, family = nb, #places notes at the deciles of sample ages
                   method = "REML")

gam.check(diatom_gam_I)
draw(diatom_gam_I)

#Compare different model fits using AIC
AIC_table <- AIC(diatom_gam_S, diatom_gam_I)%>%
  rownames_to_column(var= "Model")%>%
  mutate(data_source = rep(c("diatom_data")))%>%
  group_by(data_source)%>%
  mutate(deltaAIC = AIC - min(AIC))%>%
  ungroup()%>%
  dplyr::select(-data_source)%>%
  mutate_at(.vars = vars(df,AIC, deltaAIC), 
            .funs = funs(round,.args = list(digits=0)))

#Create synthetic data to predict over a range of ages
diat_plot_data <- with(diat_data, as_tibble(expand.grid(negAge = seq(min(diat_data$negAge), max(diat_data$negAge)),
                                                        spp = factor(levels(diat_data$spp)),
                                                        log_total_counts = mean(log_total_counts))))

# diat_modS_fit <- predict(diatom_gam_S, 
#                          newdata = diat_plot_data,
#                          se.fit = TRUE)

diat_modI_fit <- predict(diatom_gam_I,
                         newdata = diat_plot_data,
                         se.fit = TRUE)

#non-shared trends
diat_plot_data$modS_fit <- as.numeric(diat_modS_fit$fit)
diat_plot_data$modI_fit <- as.numeric(diat_modI_fit$fit)

# comparing non-shared trends
diat_plot_data <- gather(diat_plot_data, key=model, value=fit, modS_fit, modI_fit)

diat_plot_data <- mutate(diat_plot_data, se= c(as.numeric(diat_modS_fit$se.fit),
                                               as.numeric(diat_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))

# For only one model
diat_plot_data <- gather(diat_plot_data, key=model, value=fit, modI_fit)
diat_plot_data <- mutate(diat_plot_data, se= c(as.numeric(diat_modI_fit$se.fit)),
                         upper = exp(fit + (2 * se)),
                         lower = exp(fit - (2 * se)),
                         fit   = exp(fit))
diat_plot_model_labels <- paste("Model", c("I"))

#Plot the model output for non-shared trends, with means plus standard deviations for each model.
diat_plot_model_labels <- paste("Model", c("S", "I"))
diat_plot_model_labels <- factor(diat_plot_model_labels, levels = diat_plot_model_labels)

#non-shared trends
theme_set(theme_bw())
theme_update(panel.grid = element_blank())

diat_plot <- ggplot(diat_plot_data) +
  facet_wrap(~spp, nrow = 4,scales = "free_y")+
  geom_ribbon(aes(x=negAge,
                  ymin = lower,
                  ymax = upper,
                  fill = model),
              alpha=0.2)+
  geom_point(data= diat_data, aes(x = negAge, y = count), size=0.06) +
  geom_line(aes(x = negAge, y = fit, color = model))+
  labs(y = "Absolute counts", x = "Age (cal yr BP)") +
  scale_fill_brewer(name = "", palette = "Dark2",
                    labels = diat_plot_model_labels) +
  scale_colour_brewer(name = "",
                      palette = "Dark2", labels = diat_plot_model_labels)+
  theme(legend.position = "top",
        strip.text = element_text(size=10))

diat_plot

## save model results for later use
write.csv(diat_plot_data, "outputs/diatoms-HGAMs-fitted-values.csv", row.names = FALSE)

```

##Derivatives and posterior distribution simulation
```{r diatom derivatives}
# Eric's code
#Function for calculating first derivatives of time series given a before,
#after, and delta step size
calc_1st_deriv = function(fit_before, fit_after,delta) {
  (fit_after-fit_before)/(2*delta)
}

set.seed(10) #set a seed so this is repeatable
n_sims = 250

n_length = 200

years <- seq(min(diat_plot_data$negAge),
             max(diat_plot_data$negAge),
             length.out = n_length)

diff(years)
#model <- diatom_gam_S

model <- diatom_gam_I

#pred <- diat_modS_fit
pred <- diat_modI_fit

# Generate multivariate normal simulations of model coefficients
random_coefs <- t(rmvn(n_sims, mu = coef(model),V = vcov(model)))

confint_sims <- crossing(spp=unique(diat_plot_data$spp),
                         negAge = seq(min(diat_plot_data$negAge),
                                      max(diat_plot_data$negAge),
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
write.csv(deriv_summaries, "outputs/diatom-derivatives.csv", row.names = FALSE)
