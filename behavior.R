# this script creates figures of the behavior results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on behavior metrics, 2. effects of larval duration on behavior metrics) across larval period

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)

setwd("~/Desktop/R Working Directory")
# read in Rdata output that has compiled all data from EthoVision outputs
load("EthoVision_CompileData.RData")

# COMPILE DATASETS: Prepare Datasets  -----------------------
setwd("~/Desktop/R Working Directory/Databases")

#read in metamorphosis timing log database
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") 

# subset databases to only include Rana sylvatica
ev.data.summary.30 = ev.data.summary.30[ev.data.summary.30$genus.species.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS",]

# change "control" to "low density"
ev.data.summary.30$treatment = as.character(ev.data.summary.30$treatment)
ev.data.summary.30$treatment[ev.data.summary.30$treatment == "control"] = "low density"
ev.data.summary.30$treatment = factor(ev.data.summary.30$treatment)

pres.abs.summ.openfield$treatment = as.character(pres.abs.summ.openfield$treatment)
pres.abs.summ.openfield$treatment[pres.abs.summ.openfield$treatment == "control"] = "low density"
pres.abs.summ.openfield$treatment = factor(pres.abs.summ.openfield$treatment)

devo.data$treatment[devo.data$treatment == "control"] = "low density"

# change column classes
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)

pres.abs.summ.openfield$clutch = factor(pres.abs.summ.openfield$clutch)
pres.abs.summ.openfield$larv.tank.id = factor(pres.abs.summ.openfield$larv.tank.id)

# subset devo dataset to only include first six individuals
devo.data = devo.data[devo.data$first.six == "yes",]

#temporarily removed NAs (need to figure out why they exist)
ev.data.summary.30 = ev.data.summary.30[is.na(ev.data.summary.30$stim.time.cat.30) == FALSE,]

# create unique tank id for ev.data.summary.30, and devo.data
ev.data.summary.30$individual.id = as.character(ev.data.summary.30$individual.id)
ev.data.summary.30$larv.tank.id = NA
for(i in 1:nrow(ev.data.summary.30)){
  ev.data.summary.30$larv.tank.id[i] = strsplit(ev.data.summary.30$individual.id[i], split = "_")[[1]][3]
}
ev.data.summary.30$unique.id = paste(ev.data.summary.30$genus.species.code, ev.data.summary.30$clutch, ev.data.summary.30$larv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")



# ANALYZE DATA: Effect of rearing density on % frames in motion across three sampling points during larval period -----------------------

# model definition - using glmer with binomial distribution to account for proportion data

glmer.full <- glmer(movement/frames ~ treatment*week + (1|clutch:larv.tank.id), data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

glmer.nointx <- glmer(movement/frames ~ treatment + week + (1|clutch:larv.tank.id), data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

glmer.noweeks <- glmer(movement/frames ~ treatment + (1|clutch:larv.tank.id), data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

glmer.notreat <- glmer(movement/frames ~ week + (1|clutch:larv.tank.id), data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

glmer.null<- glmer(movement/frames ~ (1|clutch:larv.tank.id), data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

glm.full<- glm(movement/frames ~ treatment*week, data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], na.action = na.omit, family = binomial, weights = frames)

# model selection using likelihood ratio test
anova(glmer.full, glmer.nointx, glmer.notreat, glmer.noweeks, glmer.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher |AIC| provides better fit to the data. 
anova(glmer.full, glm.full, test = "Chisq")

# Check Model Assumptions
check_model(glmer.full)

# Final Model
Anova(glmer.full, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmer.full)
as.data.frame(exp(fixef(glmer.full))) #exponentiate the coefficient because log-transformed the response variable of mass. This gives the multiplicative factor for every one-unit increase in the independent variable.
ranef(glmer.full)
confint(glmer.full)


# ANALYZE DATA: Effect of rearing density on % arena visited across three sampling points during larval period -----------------------

temp <- pres.abs.summ.openfield %>%
  group_by(clutch, larv.tank.id, treatment, week) %>%
  summarise(n = n())

# model definition - using glmer with binomial distribution - data cannot support random effect structure of larval id nested within clutch because didn't test each larval tank at all sampling points
# OR NEED TO USE WEIGHTS = NUMBER OF VISITED SQUARES/NUMBER OF TOTAL SQUARES
# USE LET.PRESAB. POINTS TO GET AT EACH OF THESE FOR A COLUMN SO CAN USE IN THESE ANALYSES

glmer.full <- glmer(cells.presence/cells.total ~ treatment*week + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.nointx <- glmer(cells.presence/cells.total ~ treatment + week + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.noweeks <- glmer(cells.presence/cells.total ~ treatment + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.notreat <- glmer(cells.presence/cells.total ~ week + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.null<- glmer(cells.presence/cells.total ~ (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glm.full<- glm(cells.presence/cells.total ~ treatment*week, data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

# model selection using likelihood ratio test
anova(glmer.full, glmer.nointx, glmer.notreat, glmer.noweeks, glmer.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher |AIC| provides better fit to the data. 
anova(glmer.full, glm.full, test = "Chisq")

# Check Model Assumptions
check_model(glmer.full)

# Final Model
Anova(glmer.full, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmer.full)
as.data.frame(exp(fixef(glmer.full))) #exponentiate the coefficient because log-transformed the response variable of mass. This gives the multiplicative factor for every one-unit increase in the independent variable.
ranef(glmer.full)
confint(glmer.full)



# PLOT DATASETS: Effect of rearing density on % frames in motion across three sampling points during larval period -----------------------

# option 1 = plotting only open field start
plot.behav.1 <- ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=100*(movement/frames), x = week, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,100)) +
  scale_x_discrete(name = "age (weeks)")


# PLOT DATASETS: Effect of rearing density on % arena visited across three sampling points during larval period -----------------------

plot.behav.2 <- ggplot(data = pres.abs.summ.openfield, aes(y=100*prop.arena.visited, x = week, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena visited", limits = c(0,100)) +
  scale_x_discrete(name = "age (weeks)")

# PLOT DATASETS: Combine plots -----------------------
ggarrange(plot.behav.1, plot.behav.2, 
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))

