# this script creates figures of the morphology results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on morphology metrics, 2. effects of larval duration on morphology metrics) across three life stages (1. larval period, 2. at metamorphosis, 3. up to 3 months post-metamorphosis)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)
library(performance)
library(devtools)

setwd("~/Desktop/R Working Directory/Databases")


# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for metamorph morphometrics, juvenile morphometrics, and developmental timing data
morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv = read.csv("Database_Morphometrics - Froglet_Toadlet Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica and remove NAs
morph.data.mm = morph.data.mm[morph.data.mm$gs.code == "RS",]
morph.data.juv = morph.data.juv[morph.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS" & is.na(devo.data$gs.code) == FALSE,]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
morph.data.mm$treatment[morph.data.mm$treatment == "control"] = "low density"
morph.data.juv$treatment[morph.data.juv$treatment == "control"] = "low density"
devo.data$treatment[devo.data$treatment == "control"] = "low density"

#change column classes
morph.data.mm$first.six = factor(morph.data.mm$first.six, levels = c("yes", "no"))
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))

morph.data.mm$treatment = factor(morph.data.mm$treatment)
morph.data.juv$treatment = factor(morph.data.juv$treatment)
devo.data$treatment = factor(devo.data$treatment)

morph.data.mm$larv.tank.id = factor(morph.data.mm$larv.tank.id)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)

morph.data.mm$juv.tank.id = factor(morph.data.mm$juv.tank.id)
morph.data.juv$juv.tank.id = factor(morph.data.juv$juv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)

morph.data.mm$post.mm.weeks = factor(morph.data.mm$post.mm.weeks, levels = c("0"))
morph.data.juv$post.mm.weeks = factor(morph.data.juv$post.mm.weeks, levels = c("1-2", "4-6", "5-7", "8-10", "11-12", "12-14", "16-18"))

# create unique tank id for all datasets
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")

# need to change juvenile tanks with "extra" in their name to be just the name of the tank
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J002extra"] = "J002"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J003extra"] = "J003"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J005extra"] = "J005"
morph.data.juv$unique.id = paste(morph.data.juv$gs.code, morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.juv = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")

# create unique individual id for devo.data and morph.data.mm, since these are the two we have individual-level data for
morph.data.mm$unique.id.indiv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, morph.data.mm$animal.id, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, devo.data$animal.id, sep = "_")

# remove all but first six individuals from devo.data so that only those individuals are used to compute the mean
devo.data = devo.data[devo.data$first.six == "yes" & is.na(devo.data$first.six)==FALSE,]

# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.juv -----------------------

# remove overflow individuals
morph.data.juv = morph.data.juv[morph.data.juv$treatment != "overflow",]

# create column to store the developmental data
morph.data.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.mm -----------------------

# remove overflow individuals and all but first six individuals
morph.data.mm = morph.data.mm[morph.data.mm$treatment != "overflow",]
morph.data.mm = morph.data.mm[morph.data.mm$first.six == "yes",]

# create column to store the developmental data, mean.days.forelimb represents days.forelimb for each individual (rather than an average, since we actually have individual level data for this) but need to keep column name the same as in the juvenile dataset so we can combine them in the next step
morph.data.mm$mean.days.forelimb = NA

# fill developmental data column
for(i in 1:length(unique(morph.data.mm$unique.id.indiv))){
  morph.data.mm$mean.days.forelimb[morph.data.mm$unique.id.indiv == unique(morph.data.mm$unique.id.indiv)[i]] <- devo.data$days.forelimb[devo.data$unique.id.indiv ==  unique(morph.data.mm$unique.id.indiv)[i]]
}


# COMPILE DATASETS: Remove double-sampled individuals morph.data.juv  -----------------------
# some weeks had an individual measured twice - for example, if the individual was measured for both morphometrics and ephys, but only one of these should be considered since it's a repeat measure on an (unknown) individual

# only graph certain weeks (0, 1-2, 4-6, 8-10, 12-14) and when individual not measured twice within that week
morph.data.juv = morph.data.juv[morph.data.juv$post.mm.weeks != "5-7" & morph.data.juv$post.mm.weeks != "11-12",]
morph.data.juv = morph.data.juv[morph.data.juv$repeat.measure == "no",]

# COMPILE DATASETS: Combine morph.data.mm and morph.data.juv  -----------------------
morph.data.mm.juv <- rbind(morph.data.mm[ , c("date.measured",
                                             "age.weeks",
                                             "post.mm.weeks",
                                             "post.mm.sampling",
                                             "gs",
                                             "gs.code",
                                             "clutch",
                                             "treatment",
                                             "treatment.code",
                                             "juv.tank.id",
                                             "mass.g",
                                             "svl.mm",
                                             "r.forelimb.mm",    
                                              "r.tibia.mm",
                                             "r.thigh.mm",
                                             "notes.collection",
                                             "unique.id",
                                             "mean.days.forelimb")], 
                           morph.data.juv[, c("date.measured",
                                             "age.weeks",
                                             "post.mm.weeks",
                                             "post.mm.sampling",
                                             "gs",
                                             "gs.code",
                                             "clutch",
                                             "treatment",
                                             "treatment.code",
                                             "juv.tank.id",
                                             "mass.g",
                                             "svl.mm",
                                             "r.forelimb.mm",    
                                             "r.tibia.mm",
                                             "r.thigh.mm",
                                             "notes.collection",
                                             "unique.id",
                                             "mean.days.forelimb")])


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on mass at metamorphosis ---------------------

# model definition - using glmer with log-link function to keep mass in original units.

glmer.full <- glmer(mass.g ~ treatment + water.level.reduc + scale(mean.days.forelimb) + (1|clutch:larv.tank.id), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.notreat <- glmer(mass.g ~ water.level.reduc + scale(mean.days.forelimb) + (1|clutch:larv.tank.id), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log')) 

glmer.noreduc <- glmer(mass.g ~ treatment + scale(mean.days.forelimb) + (1|clutch:larv.tank.id), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.nodays <- glmer(mass.g ~ treatment + water.level.reduc + (1|clutch:larv.tank.id), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.null<- glmer(mass.g ~ (1|clutch:juv.tank.id), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log'))

glm.full<- glm(mass.g ~ treatment + water.level.reduc + scale(mean.days.forelimb), data = na.omit(morph.data.mm), na.action = na.omit, family = gaussian(link = 'log'))

# model selection using likelihood ratio test
anova(glmer.full, glmer.notreat, glmer.noreduc, glmer.nodays, glmer.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher |AIC| provides better fit to the data. 
anova(glmer.full, glm.full, test = "Chisq")

# Check Model Assumptions
check_model(glmer.full)

# Final Model
Anova(glmer.full, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmer.full)
as.data.frame(exp(fixef(glmer.full))) #exponentiate the coefficient because log-transformed the response variable of mass. This gives the multiplicative factor for every one-unit increase in the independent variable.
ranef(glmer.full)
confint(glmer.full)



# ANALYZE DATA: Effect of rearing density and larval duration on mass at and after metamorphosis ---------------------

# model definition - using glmer with log-link function to keep mass in original units.

glmer.full <- glmer(mass.g ~ treatment + scale(mean.days.forelimb)*post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.nointx <- glmer(mass.g ~ treatment + scale(mean.days.forelimb) + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glmer.nodays <- glmer(mass.g ~ treatment + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glmer.notreat <- glmer(mass.g ~ scale(mean.days.forelimb) + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log')) 

glmer.noweeks <- glmer(mass.g ~ treatment + scale(mean.days.forelimb) + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log')) 

glmer.null<- glmer(mass.g ~ (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glm.full<- glm(mass.g ~ treatment + scale(mean.days.forelimb)*post.mm.weeks, data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

# model selection using likelihood ratio test
anova(glmer.full, glmer.nointx, glmer.nodays, glmer.notreat, glmer.noweeks, glmer.null, test="Chisq")

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

#if interaction significant and fixed effect not...
#mean.days.forelimb doesn't have an effect on its own (unlike week) but it does have an effect for later weeks


# ANALYZE DATA: Effect of rearing density, treatment, and larval duration on snout-vent length at and after metamorphosis ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but ML does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. have to specify data = na.omit(morph.data.mm.juv) because some NAs in mass.g so model fitted without mass.g will have issues when we compare to other models
glmer.full <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb)*post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.nointx <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb) + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glmer.nodays <- glmer(svl.mm ~ treatment + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glmer.notreat <- glmer(svl.mm ~ scale(mean.days.forelimb) + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log')) 

glmer.noweeks <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb) + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log')) 

glmer.null<- glmer(svl.mm ~ (1|post.mm.weeks) + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glm.full<- glm(svl.mm ~ treatment + scale(mean.days.forelimb)*post.mm.weeks, data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))


# model selection using likelihood ratio test
anova(glmer.full, glmer.nointx, glmer.nodays, glmer.notreat, glmer.noweeks, glmer.null, test="Chisq")

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


# ANALYZE DATA: Correlation between morphology metrics at and after metamorphosis -----------------------

# create a dataframe of the correlation between each metric and mass for each rearing density and post-metamorphic age
for(i in 1:length(unique(morph.data.mm.juv$post.mm.weeks))){
  for(j in 1:length(unique(morph.data.mm.juv$treatment))){
    temp1 = morph.data.mm.juv[morph.data.mm.juv$post.mm.weeks == unique(morph.data.mm.juv$post.mm.weeks)[i] & morph.data.mm.juv$treatment == unique(morph.data.mm.juv$treatment)[j],]
    
    corr.svl <- cor.test(temp1$mass.g, temp1$svl.mm, 
                         method = "pearson", na.action = na.omit)
    corr.forelimb <- cor.test(temp1$mass.g, temp1$r.forelimb.mm, 
                              method = "pearson", na.action = na.omit)
    corr.tibia <- cor.test(temp1$mass.g, temp1$r.tibia.mm, 
                           method = "pearson", na.action = na.omit)
    corr.thigh <- cor.test(temp1$mass.g, temp1$r.thigh.mm, 
                           method = "pearson", na.action = na.omit)
    
    temp2 = data.frame(metric = c("svl.mm", "r.forelimb.mm", "r.tibia.mm", "r.thigh.mm"),
                       treatment = unique(temp1$treatment)[1],
                       post.mm.weeks = unique(morph.data.mm.juv$post.mm.weeks)[i],
                       type = c("svl", "forelimb", "tibia", "thigh"),
                       t = c(corr.svl$statistic, corr.forelimb$statistic, corr.tibia$statistic, corr.thigh$statistic),
                       df = c(corr.svl$parameter, corr.forelimb$parameter, corr.tibia$parameter, corr.thigh$parameter),
                       p.value = c(corr.svl$p.value, corr.forelimb$p.value, corr.tibia$p.value, corr.thigh$p.value),
                       estimate = c(corr.svl$estimate, corr.forelimb$estimate, corr.tibia$estimate, corr.thigh$estimate),
                       conf.int2.5 = c(corr.svl$conf.int[1], corr.forelimb$conf.int[1], corr.tibia$conf.int[1], corr.thigh$conf.int[1]),
                       conf.int97.5 = c(corr.svl$conf.int[2], corr.forelimb$conf.int[2], corr.tibia$conf.int[2], corr.thigh$conf.int[2])
    )
    
    if(i == 1 & j == 1){
      morph.data.corr = temp2}else{
        morph.data.corr = rbind(morph.data.corr, temp2)
      }
    rm(temp1, temp2)
  }
}

# model definition
glmer.full <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb)*post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glm.full <- glm(estimate ~ treatment + post.mm.weeks + type , data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))
glm.notreat <- glm(estimate ~ post.mm.weeks + type , data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))
glm.notype <- glm(estimate ~ treatment + post.mm.weeks, data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))
glm.null <- glm(estimate ~ post.mm.weeks, data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

# model selection using likelihood ratio test
anova(glm.full, glm.notreat, glm.notype, glm.null, test="Chisq")

# Check Model Assumptions
check_model(glm.notype)

# Final Model
summary(glm.notype)
Anova(lm.notype, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
confint(lm.notype)

# PLOT DATASETS: Effect of rearing density on morphometrics at and after metamorphosis for first six individuals from each tank -----------------------
# option 1 = x-y plot with summarized mean and +/- 1 se for all metrics

# create vector of morphometrics and y-axis labels
metrics = colnames(morph.data.juv)[24:28]
yaxis.names = c("mass (g)", "snout-vent length (mm)", "forearm length (mm)", "tibia length (mm)", "thigh length (mm)")

# create empty list to fill with morphometrics plot
plotList.morph = vector(mode = "list", length = length(metrics))

# fill plot list with each metric for mass and svl
for(i in 1:length(metrics[c(1,2)])){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics[i]]], x = post.mm.weeks, color = treatment)) + 
    #facet_grid(rows = vars(treatment)) +
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
          legend.text = element_text(size=16),
          axis.text.x=element_text(size=14, color = "black"), 
          axis.text.y=element_text(size=14, color = "black"), 
          axis.title.x=element_text(size=14, color = "black"), 
          axis.title.y = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(name = yaxis.names[i]) +
    scale_x_discrete(name = "post-metamorphic age (weeks)")
}

# fill plot list with each metric for remainder of morphometrics
for(i in 3:length(metrics)){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics[i]]], x = post.mm.weeks, color = treatment)) + 
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
          legend.text = element_text(size=14),
          axis.text.x=element_text(size=11, color = "black"), 
          axis.text.y=element_text(size=11, color = "black"), 
          axis.title.x=element_text(size=11, color = "black"), 
          axis.title.y = element_text(size=11),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(name = yaxis.names[i], limits = c(0,13)) +
    scale_x_discrete(name = "post-metamorphic age (weeks)")
}

# add the legend as the final plot within plot list so that it can be graphed within the grid
plotList.morph[[length(plotList.morph) + 1]] <- as_ggplot(get_legend(plotList.morph))

# create panel plot with only mass and svl data across all sampling points
ggarrange(plotlist = plotList.morph[c(1,2)], 
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d", "e", ""),
          font.label = list(size = 20, color = "black"))

# create panel plot with all other morphometrics data (not mass or svl) across all sampling points
ggarrange(plotlist = plotList.morph[c(3,4,5)], 
          ncol = 1,
          nrow = 3,
          common.legend = TRUE,
          legend = "right",
          labels = c("a", "b", "c", "d", "e", ""),
          font.label = list(size = 16, color = "black"))



# PLOT DATASETS: Effect of rearing density on correlation between morphology metrics at and after metamorphosis -----------------------

# create empty list to fill with morphometrics plot
plotList.morph = vector(mode = "list", length = length(metrics[-1]))

# fill plot list with each metric
for(i in 2:length(metrics)){
  plotList.morph[[i-1]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics[i]]], x = mass.g, color = post.mm.weeks)) + 
    #facet_grid(rows = vars(treatment)) +
    geom_abline(slope = 1, intercept = 0, na.rm = FALSE, show.legend = NA, linetype = 2) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
    scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
    scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text.x=element_text(size=12, color = "black"), 
          axis.text.y=element_text(size=12, color = "black"), 
          axis.title.x=element_text(size=12, color = "black"), 
          axis.title.y = element_text(size=12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    expand_limits(y = 0) +
    scale_y_continuous(name = yaxis.names[i]) +
    scale_x_continuous(name = "mass (g)")
}

# create panel plot with all morphometrics data across all sampling points
ggarrange(plotlist = plotList.morph,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))


# PLOT DATASETS: Correlation matrix between morphology metrics at and after metamorphosis -----------------------

morph.data.corr$metric = factor(morph.data.corr$metric, levels = c("svl.mm", "r.forelimb.mm", "r.tibia.mm", "r.thigh.mm"))

ggplot(data = morph.data.corr, aes(y=estimate, x = metric, color = post.mm.weeks)) + 
  facet_grid(rows = vars(treatment)) +
  geom_point(position=position_jitterdodge(), size = 5, alpha = 1) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=14, color = "black"), 
        axis.text.y=element_text(size=14, color = "black"), 
        axis.title.x=element_text(size=14, color = "black"), 
        axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "correlation estimate") +
  scale_x_discrete(name = "", labels = c("snout-vent", "forearm", "tibia", "thigh"))



# PLOT DATASETS: Effect of developmental speed on morphometrics at and after metamorphosis ---------------------
# option 1 = x-y plot with summarized mean and +/- 1 se for all metrics

# create empty list to fill with morphometrics plot
plotList.morph = vector(mode = "list", length = length(metrics))

# fill plot list with each metric
for(i in 1:length(metrics)){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics[i]]], x = mean.days.forelimb, color = post.mm.weeks, fill = post.mm.weeks)) + 
    #facet_grid(rows = vars(post.mm.weeks)) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_smooth(method="glm",
                method.args = list(family = gaussian(link = 'log')), 
                mapping = aes(color=as.factor(post.mm.weeks)), 
                se=T) +
    scale_color_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
    scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3]))+
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text.x=element_text(size=14, color = "black"), 
          axis.text.y=element_text(size=14, color = "black"), 
          axis.title.x=element_text(size=14, color = "black"), 
          axis.title.y = element_text(size=14),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    scale_y_continuous(name = yaxis.names[i], n.breaks = 5) +
    scale_x_continuous(name = "larval duration (days)", limits = c(45, 75), breaks = seq(45, 75, by = 5))
}

# add the legend as the final plot within plot list so that it can be graphed within the grid
plotList.morph[[length(plotList.morph) + 1]] <- as_ggplot(get_legend(plotList.morph))

# create panel plot with all morphometrics data across all sampling points
ggarrange(plotlist = plotList.morph, 
          legend = "none",
          labels = c("a", "b", "c", "d", "e", ""),
          font.label = list(size = 20, color = "black"))

# create panel plot with only mass and svl data across all sampling points
ggarrange(plotlist = plotList.morph[1:2],
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))


# SUMMARY TABLES: Morphological metrics at and after metamorphosis ---------------------

#reset working directory to be outputs folder
setwd("~/Desktop/R Working Directory/Outputs")

# create summary table for morphometrics that doesn't include rearing density
morph.data.summ = morph.data.mm.juv %>%
  group_by(post.mm.weeks) %>%
  summarise(mean.mass = round(mean(mass.g, na.rm = TRUE), 4),
            sd.mass = round(sd(mass.g, na.rm = TRUE), 4),
            mean.svl = round(mean(svl.mm, na.rm = TRUE), 2),
            sd.svl = round(sd(svl.mm, na.rm = TRUE), 2),
            mean.forearm = round(mean(r.forelimb.mm, na.rm = TRUE), 2),
            sd.forearm = round(sd(r.forelimb.mm, na.rm = TRUE), 2),
            mean.tibia = round(mean(r.tibia.mm, na.rm = TRUE), 2),
            sd.tibia = round(sd(r.tibia.mm, na.rm = TRUE),2),
            mean.thigh = round(mean(r.thigh.mm, na.rm = TRUE),2),
            sd.thigh = round(sd(r.thigh.mm, na.rm = TRUE),2)
            )

# create summary table for correlation matrix that doesn't include rearing density
for(i in 1:length(unique(morph.data.mm.juv$post.mm.weeks))){
    temp1 = morph.data.mm.juv[morph.data.mm.juv$post.mm.weeks == unique(morph.data.mm.juv$post.mm.weeks)[i],]
    
    corr.svl <- cor.test(temp1$mass.g, temp1$svl.mm, 
                         method = "pearson", na.action = na.omit)
    corr.forelimb <- cor.test(temp1$mass.g, temp1$r.forelimb.mm, 
                              method = "pearson", na.action = na.omit)
    corr.tibia <- cor.test(temp1$mass.g, temp1$r.tibia.mm, 
                           method = "pearson", na.action = na.omit)
    corr.thigh <- cor.test(temp1$mass.g, temp1$r.thigh.mm, 
                           method = "pearson", na.action = na.omit)
    
    temp2 = data.frame(metric = c("svl.mm", "r.forelimb.mm", "r.tibia.mm", "r.thigh.mm"),
                       post.mm.weeks = unique(morph.data.mm.juv$post.mm.weeks)[i],
                       type = c("svl", "forelimb", "tibia", "thigh"),
                       t = c(corr.svl$statistic, corr.forelimb$statistic, corr.tibia$statistic, corr.thigh$statistic),
                       df = c(corr.svl$parameter, corr.forelimb$parameter, corr.tibia$parameter, corr.thigh$parameter),
                       p.value = c(corr.svl$p.value, corr.forelimb$p.value, corr.tibia$p.value, corr.thigh$p.value),
                       estimate = round(c(corr.svl$estimate, corr.forelimb$estimate, corr.tibia$estimate, corr.thigh$estimate), 2),
                       conf.int2.5 = round(c(corr.svl$conf.int[1], corr.forelimb$conf.int[1], corr.tibia$conf.int[1], corr.thigh$conf.int[1]),2),
                       conf.int97.5 = round(c(corr.svl$conf.int[2], corr.forelimb$conf.int[2], corr.tibia$conf.int[2], corr.thigh$conf.int[2]), 2)
    )
    
    if(i == 1){
      morph.data.corr.summ = temp2}else{
        morph.data.corr.summ = rbind(morph.data.corr.summ, temp2)
      }
    rm(temp1, temp2)
}
morph.data.corr.summ$conf.int = paste(morph.data.corr.summ$conf.int2.5, morph.data.corr.summ$conf.int97.5, sep = "-")
write.table(morph.data.corr.summ, file = "MorphDataCorrSumm.txt", sep = ",")
