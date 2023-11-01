# this script creates figures of the behavior results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on behavior metrics, 2. effects of larval duration on behavior metrics) across larval period

#helpful references for interpreting odds ratios
#https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-how-do-i-interpret-odds-ratios-in-logistic-regression/
#https://www.bookdown.org/rwnahhas/RMPH/blr-interp.html

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)
library(performance)

setwd("~/Desktop/R Working Directory")
# read in Rdata output that has compiled all data from EthoVision outputs
load("EthoVision_CompileData.RData")

# COMPILE DATASETS: Prepare Datasets  -----------------------
setwd("~/Desktop/R Working Directory/Databases")

#read in databases
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") 
morph.data.tad = read.csv("Database_Morphometrics - Tadpole Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
behav.data.tad = read.csv("Database_Behavior - Startle - Tadpole.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
ev.data.summary.30 = ev.data.summary.30[ev.data.summary.30$genus.species.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS",]
morph.data.tad = morph.data.tad[morph.data.tad$gs.code == "RS",]
behav.data.tad = behav.data.tad[behav.data.tad$gs.code == "RS",]


# subset ev.data.summary.30 to only include open-field data
ev.data.summary.30 = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",]

# subset morphological database to only include behavior measurements
morph.data.tad = morph.data.tad[morph.data.tad$data.type == "behavior",]

# change "control" to "low density"
ev.data.summary.30$treatment = as.character(ev.data.summary.30$treatment)
ev.data.summary.30$treatment[ev.data.summary.30$treatment == "control"] = "low density"
ev.data.summary.30$treatment = factor(ev.data.summary.30$treatment)

pres.abs.summ.openfield$treatment = as.character(pres.abs.summ.openfield$treatment)
pres.abs.summ.openfield$treatment[pres.abs.summ.openfield$treatment == "control"] = "low density"
pres.abs.summ.openfield$treatment = factor(pres.abs.summ.openfield$treatment)

devo.data$treatment[devo.data$treatment == "control"] = "low density"

morph.data.tad$treatment[morph.data.tad$treatment == "control"] = "low density"

behav.data.tad$treatment[behav.data.tad$treatment == "control"] = "low density"

# change column classes
ev.data.summary.30$week = factor(ev.data.summary.30$week, ordered = TRUE, levels = c("3", "6", "9"))

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

# create unique tank id for ev.data.summary.30, pres.abs.summ.openfield, morph.data.tad, behav.data.tad, and devo.data
ev.data.summary.30$individual.id = as.character(ev.data.summary.30$individual.id)
ev.data.summary.30$larv.tank.id = NA
ev.data.summary.30$animal.id.num = NA
ev.data.summary.30$date.photo = NA
for(i in 1:nrow(ev.data.summary.30)){
  ev.data.summary.30$larv.tank.id[i] = strsplit(ev.data.summary.30$individual.id[i], split = "_")[[1]][3]
  ev.data.summary.30$animal.id.num[i] = as.numeric(strsplit(ev.data.summary.30$individual.id[i], split = "_")[[1]][6])
  ev.data.summary.30$date.photo[i] = strsplit(ev.data.summary.30$individual.id[i], split = "_")[[1]][4]
}
ev.data.summary.30$unique.id = paste(ev.data.summary.30$genus.species.code, ev.data.summary.30$clutch, ev.data.summary.30$larv.tank.id, ev.data.summary.30$date.photo, ev.data.summary.30$animal.id.num, sep = "_")

pres.abs.summ.openfield$individual.id = as.character(pres.abs.summ.openfield$individual.id)
pres.abs.summ.openfield$larv.tank.id = NA
pres.abs.summ.openfield$animal.id.num = NA
pres.abs.summ.openfield$date.photo = NA
for(i in 1:nrow(pres.abs.summ.openfield)){
  pres.abs.summ.openfield$larv.tank.id[i] = strsplit(pres.abs.summ.openfield$individual.id[i], split = "_")[[1]][3]
  pres.abs.summ.openfield$animal.id.num[i] = as.numeric(strsplit(pres.abs.summ.openfield$individual.id[i], split = "_")[[1]][6])
  pres.abs.summ.openfield$date.photo[i] = strsplit(pres.abs.summ.openfield$individual.id[i], split = "_")[[1]][4]
}
pres.abs.summ.openfield$unique.id = paste(pres.abs.summ.openfield$genus.species.code, pres.abs.summ.openfield$clutch, pres.abs.summ.openfield$larv.tank.id, pres.abs.summ.openfield$date.photo, pres.abs.summ.openfield$animal.id.num, sep = "_")


devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")

morph.data.tad$unique.id = paste(morph.data.tad$gs.code, morph.data.tad$clutch, morph.data.tad$larv.tank.id, morph.data.tad$date.photo , morph.data.tad$animal.id.num , sep = "_")

behav.data.tad$unique.id = paste(behav.data.tad$gs.code, behav.data.tad$clutch, behav.data.tad$larv.tank.id, behav.data.tad$video.date, behav.data.tad$animal.id.num , sep = "_")


# COMPILE DATASETS: Combine morphological and gosner stage data with ev.data.summary.30 and pres.abs.summ.openfield  -----------------------
ev.data.summary.30$tl.cm = NA
ev.data.summary.30$bl.cm = NA
ev.data.summary.30$hw.cm = NA
ev.data.summary.30$bw.cm = NA
ev.data.summary.30$tmw.cm = NA
ev.data.summary.30$avg.gosner.stage = NA

for(i in 1:length(unique(ev.data.summary.30$unique.id))){
  ev.data.summary.30$tl.cm[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = morph.data.tad$tl.cm[morph.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
  
  ev.data.summary.30$bl.cm[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = morph.data.tad$bl.cm[morph.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
  
  ev.data.summary.30$hw.cm[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = morph.data.tad$hw.cm[morph.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
  
  ev.data.summary.30$bw.cm[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = morph.data.tad$bw.cm[morph.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
  
  ev.data.summary.30$tmw.cm[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = morph.data.tad$tmw.cm[morph.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
  
  ev.data.summary.30$avg.gosner.stage[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] = behav.data.tad$avg.gosner.stage[behav.data.tad$unique.id == unique(ev.data.summary.30$unique.id)[i]] 
}


pres.abs.summ.openfield$tl.cm = NA
pres.abs.summ.openfield$bl.cm = NA
pres.abs.summ.openfield$hw.cm = NA
pres.abs.summ.openfield$bw.cm = NA
pres.abs.summ.openfield$tmw.cm = NA
pres.abs.summ.openfield$avg.gosner.stage = NA

for(i in 1:length(unique(pres.abs.summ.openfield$unique.id))){
  pres.abs.summ.openfield$tl.cm[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = morph.data.tad$tl.cm[morph.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] 
  
  pres.abs.summ.openfield$bl.cm[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = morph.data.tad$bl.cm[morph.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] 
  
  pres.abs.summ.openfield$hw.cm[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = morph.data.tad$hw.cm[morph.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]]
  
  pres.abs.summ.openfield$bw.cm[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = morph.data.tad$bw.cm[morph.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]]
  
  pres.abs.summ.openfield$tmw.cm[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = morph.data.tad$tmw.cm[morph.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]]
  
  pres.abs.summ.openfield$avg.gosner.stage[pres.abs.summ.openfield$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] = behav.data.tad$avg.gosner.stage[behav.data.tad$unique.id == unique(pres.abs.summ.openfield$unique.id)[i]] 
}


# ANALYZE DATA: Effect of rearing density on % frames in motion across three sampling points during larval period -----------------------
#test correlations of fixed effects
cor.test(ev.data.summary.30$week, ev.data.summary.30$bl.cm,
         method = "pearson", na.action = na.omit)
cor.test(ev.data.summary.30$week, ev.data.summary.30$tl.cm,
         method = "pearson", na.action = na.omit)
cor.test(ev.data.summary.30$week, ev.data.summary.30$avg.gosner.stage,
         method = "pearson", na.action = na.omit)

cor.test(ev.data.summary.30$bl.cm, ev.data.summary.30$avg.gosner.stage,
         method = "pearson", na.action = na.omit)
cor.test(ev.data.summary.30$tl.cm, ev.data.summary.30$avg.gosner.stage,
         method = "pearson", na.action = na.omit)

cor.test(ev.data.summary.30$bl.cm, ev.data.summary.30$tl.cm,
         method = "pearson", na.action = na.omit)

# model definition - using glmer with binomial distribution to account for proportion data
glmer.fullintxn <- glmer(movement/frames ~ treatment*week + treatment*bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.weekintxn <- glmer(movement/frames ~ treatment*week + bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.sizeintxn <- glmer(movement/frames ~ treatment + week + treatment*bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.nointx <- glmer(movement/frames ~ treatment + week + bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.noweeks <- glmer(movement/frames ~ treatment + bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.notreat <- glmer(movement/frames ~ week + bl.cm + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.nosize <- glmer(movement/frames ~ treatment + week + (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glmer.null<- glmer(movement/frames ~ (1|clutch:larv.tank.id), data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

glm.full<- glm(movement/frames ~ treatment*week + bl.cm, data = ev.data.summary.30, na.action = na.omit, family = binomial, weights = frames)

# model selection using likelihood ratio test
anova(glmer.fullintxn, glmer.weekintxn, glmer.sizeintxn, glmer.nointx, glmer.notreat, glmer.nosize, glmer.noweeks, glmer.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher |AIC| provides better fit to the data. 
anova(glmer.fullintxn, glm.full, test = "Chisq")

# Check Model Assumptions
check_model(glmer.fullintxn)

# Final Model
Anova(glmer.fullintxn, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmer.fullintxn)
as.data.frame(exp(fixef(glmer.fullintxn))) #exponentiate the coefficient to transform log odds to an odds ratio - good explanation here: https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-how-do-i-interpret-odds-ratios-in-logistic-regression/
ranef(glmer.full)
confint(glmer.full)

#for high density individuals, the odds ratio is XX for a one-unit increase in week
exp(fixef(glmer.fullintxn)[3])
#for low density individuals, the odds ratio is XX for a one-unit increase in week
exp(fixef(glmer.fullintxn)[3] + fixef(glmer.fullintxn)[6])

#for high density individuals, the odds ratio is XX for a one-unit increase in body length
exp(fixef(glmer.fullintxn)[5])
#for low density individuals, the odds ratio is XX for a one-unit increase in week
exp(fixef(glmer.fullintxn)[5] + fixef(glmer.fullintxn)[8])


# ANALYZE DATA: Effect of rearing density on % arena visited across three sampling points during larval period -----------------------

temp <- pres.abs.summ.openfield %>%
  group_by(clutch, larv.tank.id, treatment, week) %>%
  summarise(n = n())

# model definition - using glmer with binomial distribution - data cannot support random effect structure of larval id nested within clutch because didn't test each larval tank at all sampling points
# OR NEED TO USE WEIGHTS = NUMBER OF VISITED SQUARES/NUMBER OF TOTAL SQUARES
# USE LET.PRESAB. POINTS TO GET AT EACH OF THESE FOR A COLUMN SO CAN USE IN THESE ANALYSES
pres.abs.summ.openfield$week = as.numeric(pres.abs.summ.openfield$week)

glmer.fullintxn <- glmer(cells.presence/cells.total ~ treatment*scale(week) + treatment*bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.weekintxn <- glmer(cells.presence/cells.total ~ treatment*scale(week) + bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.sizeintxn <- glmer(cells.presence/cells.total ~ treatment + scale(week) + treatment*bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.nointx <- glmer(cells.presence/cells.total ~ treatment + week + bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.noweeks <- glmer(cells.presence/cells.total ~ treatment + bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.nosize <- glmer(cells.presence/cells.total ~ treatment + week + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.notreat <- glmer(cells.presence/cells.total ~ week + bl.cm + (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glmer.null<- glmer(cells.presence/cells.total ~ (1|clutch:larv.tank.id), data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

glm.full<- glm(cells.presence/cells.total ~ treatment*week + treatment*bl.cm, data = pres.abs.summ.openfield, na.action = na.omit, family = binomial, weights = cells.total)

# model selection using likelihood ratio test
anova(glmer.fullintxn, glmer.weekintxn, glmer.sizeintxn, glmer.nointx, glmer.notreat, glmer.noweeks, glmer.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher |AIC| provides better fit to the data. 
anova(glmer.fullintxn, glm.full, test = "Chisq")

# Check Model Assumptions
check_model(glmer.full)

# Final Model
Anova(glmer.fullintxn, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmer.fullintxn)
as.data.frame(exp(fixef(glmer.full))) #exponentiate the coefficient because log-transformed the response variable of mass. This gives the multiplicative factor for every one-unit increase in the independent variable.
ranef(glmer.full)
confint(glmer.full)



# PLOT DATASETS: Effect of rearing density on % frames in motion across three sampling points during larval period -----------------------

# option 1 = plotting only open field start
plot.behav.1 <- ggplot(data = ev.data.summary.30, aes(y=100*(movement/frames), x = week, color = treatment)) + 
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
  scale_x_continuous(name = "age (weeks)", breaks = c(3,6,9))

ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=100*(movement/frames), x = avg.gosner.stage, color = treatment)) + 
  facet_grid(rows= vars(week)) +
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
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
  scale_x_continuous(name = "average gosner stage")

ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=100*(movement/frames), x = bl.cm, color = treatment)) + 
  facet_grid(rows= vars(week)) +
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
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
  scale_x_continuous(name = "body length (cm)", limits = c(0.5, 1.5))

ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=avg.gosner.stage, x = as.numeric(as.character(week)), color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
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
  scale_y_continuous(name = "average gosner stage") +
  scale_x_continuous(name = "pre-metamorphic age (weeks)", limits = c(3,9), breaks = c(3,6,9))


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


# PLOT DATASETS: Difference in body size (cm) across three sampling points during larval period -----------------------

#correlation test
cor.test(ev.data.summary.30$bl.cm, ev.data.summary.30$tl.cm,
         method = "pearson", na.action = na.omit)

# total length
plot.size.tl <- ggplot(data = ev.data.summary.30, aes(y=tl.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=tl.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "total length (cm)") +
  scale_x_discrete(name = "week")

# body length
plot.size.bl <- ggplot(data = ev.data.summary.30, aes(y=bl.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=bl.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body length (cm)") +
  scale_x_discrete(name = "week")


# body length/total length
plot.size.bltl <- ggplot(data = ev.data.summary.30, aes(y=bl.cm/tl.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=bl.cm/tl.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body length/total length (cm)") +
  scale_x_discrete(name = "week")

# body width
plot.size.bw <- ggplot(data = ev.data.summary.30, aes(y=bw.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=bw.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body width (cm)") +
  scale_x_discrete(name = "week")

# head width
plot.size.hw <- ggplot(data = ev.data.summary.30, aes(y=hw.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=hw.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "head width (cm)") +
  scale_x_discrete(name = "week")


# tail-muscle width
plot.size.tmw <- ggplot(data = ev.data.summary.30, aes(y=tmw.cm, x = as.factor(week), color = treatment), show.legend = FALSE) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(aes(y=tmw.cm, x = as.factor(week))) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "tail-muscle width (cm)") +
  scale_x_discrete(name = "week")


# PLOT DATASETS: Effect of body size (cm) on % frames in motion across three sampling points during larval period -----------------------

plot.move.tl <- ggplot(data = ev.data.summary.30, aes(y=(movement/frames), x = tl.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial),
              mapping = aes(color=treatment, fill = treatment, weight = frames), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "total length (cm)")


plot.move.bl <- ggplot(data = ev.data.summary.30, aes(y=movement/frames, x = bl.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial), 
              mapping = aes(color=treatment, fill = treatment, weight = frames), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "body length (cm)")


plot.move.bw <- ggplot(data = ev.data.summary.30, aes(y=movement/frames, x = bw.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial), 
              mapping = aes(color=treatment, fill = treatment, weight = frames), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "body width (cm)")


plot.move.tmw <- ggplot(data = ev.data.summary.30, aes(y=movement/frames, x = tmw.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial), 
              mapping = aes(color=treatment, fill = treatment, weight = frames), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "tail-muscle width (cm)", limits = c(0,0.5))



# PLOT DATASETS: Effect of body size (cm) on exploration across three sampling points during larval period -----------------------

plot.explo.tl <- ggplot(data = pres.abs.summ.openfield, aes(y=(cells.presence/cells.total), x = tl.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial),
              mapping = aes(color=treatment, fill = treatment, weight = cells.total), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena explored", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "total length (cm)")


plot.explo.bl <- ggplot(data = pres.abs.summ.openfield, aes(y=(cells.presence/cells.total), x = bl.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial),
              mapping = aes(color=treatment, fill = treatment, weight = cells.total), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena explored", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "body length (cm)")


plot.explo.bw <- ggplot(data = pres.abs.summ.openfield, aes(y=(cells.presence/cells.total), x = bw.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial),
              mapping = aes(color=treatment, fill = treatment, weight = cells.total), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena explored", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "body width (cm)")


plot.explo.tmw <- ggplot(data = pres.abs.summ.openfield, aes(y=(cells.presence/cells.total), x = tmw.cm, color = treatment, show.legend = FALSE)) + 
  facet_grid(cols = vars(week)) +
  geom_point(size = 2.5, alpha = 0.7, show.legend = TRUE) +
  stat_smooth(method="glm", #fit logistic model with confidence interval
              method.args = list(family = binomial),
              mapping = aes(color=treatment, fill = treatment, weight = cells.total), 
              se=T) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x=element_text(size=14, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena explored", limits = c(0,1), labels = c(0, 25, 50, 75, 100)) +
  scale_x_continuous(name = "tail-muscle width (cm)")



# PLOT DATASETS: Combine plots -----------------------
ggarrange(plot.size.tl, plot.size.bl, plot.size.bltl, plot.size.hw, plot.size.bw, plot.size.tmw,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d", "e", "f"),
          font.label = list(size = 20, color = "black"))

#morphometrics
ggarrange(plot.move.tl, plot.move.bl, plot.move.bw, plot.move.tmw,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))

ggarrange(plot.explo.tl, plot.explo.bl, plot.explo.bw, plot.explo.tmw,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))

