# this script creates figures of the survival results for the 2023 ranid developmental paper
# includes effects of rearing density on survival across three stages (1. larval duration, at metamorphosis, after metamorphosis)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(coxme) #don't need this package if not running Cox proportional hazards models
library(car)
library(lme4)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for survival data and development data
survi.data.tad = read.csv("Database_Survivorship - Tadpole Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.exp = read.csv("Database_Survivorship - Tadpole Experiment Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.juv = read.csv("Database_Survivorship - Froglet_Toadlet Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database

# subset databases to only include Rana sylvatica
survi.data.tad = survi.data.tad[survi.data.tad$gs.code == "RS",]
survi.data.exp = survi.data.exp[survi.data.exp$gs.code == "RS",]
survi.data.juv = survi.data.juv[survi.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS" & is.na(devo.data$gs.code) == FALSE,]

# subset databases to only include weeks up to 18 (currently at 17 but need to change after October 23)
survi.data.juv = survi.data.juv[survi.data.juv$postmm.week <= 18 & is.na(survi.data.juv$postmm.week) == FALSE,]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
survi.data.tad$treatment[survi.data.tad$treatment == "control"] = "low density"
survi.data.exp$treatment[survi.data.exp$treatment == "control"] = "low density"
survi.data.juv$treatment[survi.data.juv$treatment == "control"] = "low density"
devo.data$treatment[devo.data$treatment == "control"] = "low density"

#change column classes
survi.data.tad$treatment = factor(survi.data.tad$treatment)
survi.data.tad$larv.tank.id = factor(survi.data.tad$larv.tank.id)

survi.data.exp$treatment = factor(survi.data.exp$treatment)
survi.data.exp$larv.tank.id = factor(survi.data.exp$larv.tank.id)

survi.data.juv$treatment = factor(survi.data.juv$treatment)
survi.data.juv$juv.tank.id = factor(survi.data.juv$juv.tank.id)

devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)
devo.data$days.forelimb = as.integer(devo.data$days.forelimb)
devo.data$days.forelimb.tail = as.integer(devo.data$days.forelimb.tail)

# subset devo dataset to only include first six individuals
devo.data = devo.data[devo.data$first.six == "yes",]

# create unique tank id for survi.data.tad, survi.data.exp, survi.data.juv, and devo.data
survi.data.tad$unique.id = paste(survi.data.tad$gs.code, survi.data.tad$clutch, survi.data.tad$larv.tank.id, sep = "_")
survi.data.exp$unique.id = paste(survi.data.exp$gs.code, survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_")
survi.data.juv$unique.id = paste(survi.data.juv$gs.code, survi.data.juv$clutch, survi.data.juv$juv.tank.id, sep = "_")
survi.data.juv$unique.id.wk = paste(survi.data.juv$gs.code, survi.data.juv$clutch, survi.data.juv$juv.tank.id, survi.data.juv$postmm.week, sep = "_")
devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.juv = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")


# COMPILE DATASETS: Combine weekly tadpole survivorship with experiment tadpole survivorship -----------
colnames(survi.data.exp)[11] = "photo.num"
colnames(survi.data.exp)[12] = "leth.samp.num.cumul"

survi.data.tad.exp = rbind(survi.data.tad[ , c("week",
                                               "gs",
                                               "gs.code",
                                               "clutch",
                                               "treatment",
                                               "treatment.code",
                                               "larv.tank.id",
                                               "seed.num",
                                               "photo.num",
                                               "leth.samp.num.cumul",
                                               "metamorph.num.cumul", 
                                               "prop.survi",
                                               "prop.seed.forelimb",
                                               "prop.surv.forelimb",
                                               "water.level.reduc",
                                               "unique.id")], 
                           survi.data.exp[ , c("week",
                                               "gs",
                                               "gs.code",
                                               "clutch",
                                               "treatment",
                                               "treatment.code",
                                               "larv.tank.id",
                                               "seed.num",
                                               "photo.num",
                                               "leth.samp.num.cumul",
                                               "metamorph.num.cumul", 
                                               "prop.survi",
                                               "prop.seed.forelimb",
                                               "prop.surv.forelimb",
                                               "water.level.reduc",
                                               "unique.id")])

# create unique tank id for newly created survi.data.tad.exp
survi.data.tad.exp$unique.id = paste(survi.data.tad.exp$gs.code, survi.data.tad.exp$clutch, survi.data.tad.exp$larv.tank.id, sep = "_")


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to survi.data.tad.exp -----------------------

# create column to store the developmental data, mean.days.forelimb represents days.forelimb for each individual (rather than an average, since we actually have individual level data for this) but need to keep column name the same as in the juvenile dataset so we can combine them in the next step
survi.data.tad.exp$mean.days.forelimb = NA

# fill developmental data column for survi.data.tad.exp
for(i in 1:length(unique(survi.data.tad.exp$unique.id))){
  survi.data.tad.exp$mean.days.forelimb[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(survi.data.tad.exp$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to survi.data.juv -----------------------

# create column to store the developmental data, mean.days.forelimb represents days.forelimb for each individual (rather than an average, since we actually have individual level data for this) but need to keep column name the same as in the juvenile dataset so we can combine them in the next step
survi.data.juv$mean.days.forelimb = NA

# fill developmental data column for survi.data.juv
for(i in 1:length(unique(survi.data.juv$unique.id))){
  survi.data.juv$mean.days.forelimb[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(survi.data.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical developmental data (i.e. early, mid, late developers)
survi.data.juv$devo.cat = NA

for(i in 1:length(unique(survi.data.juv$unique.id))){
  
  if(survi.data.juv$mean.days.forelimb[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] <= mean(survi.data.juv$mean.days.forelimb, na.rm = TRUE) - sd(survi.data.juv$mean.days.forelimb, na.rm = TRUE)){
    survi.data.juv$devo.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "early"
  }
  
  if(survi.data.juv$mean.days.forelimb[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] >= mean(survi.data.juv$mean.days.forelimb, na.rm = TRUE) + sd(survi.data.juv$mean.days.forelimb, na.rm = TRUE)){
    survi.data.juv$devo.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "late"
  }
}
survi.data.juv$devo.cat[is.na(survi.data.juv$devo.cat) == TRUE] = "mid"


# COMPILE DATA: Generate survivorship dataframe in terms of 0 (dead) and 1 (alive) for tadpole weeks up to end of experiment ---------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[1],]
temp3 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                   gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                   gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                   clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                   status = NA
)
temp3$status[0:(temp$photo.num + temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                     gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                     status = NA
  )
  temp2$status[1:(temp$photo.num + temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.tad.exp.longform = temp3
rm(temp3)


# COMPILE DATA: Generate survivorship dataframe in terms of 0 (dead) and 1 (alive) for juvenile weeks ---------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.juv[survi.data.juv$unique.id.wk == unique(survi.data.juv$unique.id.wk)[1],]
temp3 = data.frame(postmm.week = rep(temp$postmm.week, temp$seed.num - temp$leth.samp.num),
                   gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                   gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                   clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                   treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                   treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                   juv.tank.id = rep(temp$juv.tank.id, temp$seed.num - temp$leth.samp.num),
                   mean.days.forelimb = rep(temp$mean.days.forelimb, temp$seed.num - temp$leth.samp.num),
                   status = NA
)
temp3$status[0:(temp$live.num)] = 1
temp3$status[is.na(temp$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.juv$unique.id.wk))){
  temp = survi.data.juv[survi.data.juv$unique.id.wk == unique(survi.data.juv$unique.id.wk)[i],]
  temp2 = data.frame(postmm.week = rep(temp$postmm.week, temp$seed.num - temp$leth.samp.num),
                     gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                     juv.tank.id = rep(temp$juv.tank.id, temp$seed.num - temp$leth.samp.num),
                     mean.days.forelimb = rep(temp$mean.days.forelimb, temp$seed.num - temp$leth.samp.num),
                     status = NA
  )
  temp2$status[0:(temp$live.num)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.juv.longform = temp3
rm(temp3)


# ANALYZE DATA: Effect of water level reduction (yes/no) on larval duration ---------------------
# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but ML does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.

# model definition - setting one random effect as larval tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
glmm.full <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.notreat <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ week + water.level.reduc + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.noweek <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + water.level.reduc + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.noreduc <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.null <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glm.full <- glm(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc, data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.nointxn, glmm.notreat, glmm.noweek, glmm.noreduc, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.nointxn, glm.full, test="Chisq") 

# Final Model
Anova(glmm.nointxn, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmm.nointxn)
exp(fixef(glmm.nointxn)) #returns odds ratis instead of logit scale
exp(confint(glmm.nointxn)) #returns odds ratis instead of logit scale
ranef(glmm.nointxn)
VarCorr(glmm.nointxn)


# ANALYZE DATA: Cox Proportional Hazard Model - Effect of rearing density, week, and larv tank id nested within clutch on likelihood of survival weekly and by close-out of tanks ---------------------

# model definition - setting one random effect as week and one random effect as larval tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
cox.model.mixed.full <- coxme(Surv(status) ~ treatment + (1|week) + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model.mixed.noweek <- coxme(Surv(status) ~ treatment + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model.mixed.notank <- coxme(Surv(status) ~ treatment + (1|week), data=survi.data.tad.exp.longform)
cox.model.mixed.null <- coxme(Surv(status) ~ (1|week) + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model <- coxph(Surv(status) ~ treatment, data=survi.data.tad.exp.longform)

# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model.mixed.noweek, cox.model.mixed.notank, cox.model.mixed.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model, test="Chisq") 

# Final Model
Anova(cox.model.mixed.full, type = "II")

summary(cox.model.mixed.full)
exp(confint(cox.model.mixed.full))
VarCorr(cox.model.mixed.full)



# ANALYZE DATA: Cox Proportional Hazard Model - Effect of rearing density, larval duration, week, and juv tank id nested within clutch on likelihood of survival weekly ---------------------

# model definition - setting one random effect as week and one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16
cox.model.mixed.full <- coxme(Surv(status) ~ treatment*mean.days.forelimb + (1|postmm.week) + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.nointxn <- coxme(Surv(status) ~ treatment + mean.days.forelimb + (1|postmm.week) + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.notreat <- coxme(Surv(status) ~ mean.days.forelimb + (1|postmm.week) + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.nodays <- coxme(Surv(status) ~ treatment + (1|postmm.week) + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.noweek <- coxme(Surv(status) ~ treatment + mean.days.forelimb + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.notank <- coxme(Surv(status) ~ treatment + mean.days.forelimb + (1|postmm.week), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model.mixed.null <- coxme(Surv(status) ~ (1|postmm.week) + (1|clutch/juv.tank.id), data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])
cox.model <- coxph(Surv(status) ~ treatment + mean.days.forelimb, data=survi.data.juv.longform[survi.data.juv.longform$postmm.week < 13 & survi.data.juv.longform$postmm.week > 0,])

# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model.mixed.nointxn, cox.model.mixed.notreat, cox.model.mixed.nodays, cox.model.mixed.noweek, cox.model.mixed.notank, cox.model.mixed.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model, test="Chisq") 

# Final Model
Anova(cox.model.mixed.notreat, type = "II")
summary(cox.model.mixed.notreat)
confint(cox.model.mixed.notreat)
VarCorr(cox.model.mixed.notreat)


# ANALYZE DATA: Effect of rearing density, larval duration, week, and larv tank id nested within clutch on % surviving individuals weekly up until metamorphosis ---------------------

# model definition - setting one random effect as larval tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
glmm.full <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment*week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment*week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.notreat <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.noweek <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.null <- glmer(photo.num/(seed.num-leth.samp.num.cumul) ~ (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glm.full <- glm(photo.num/(seed.num-leth.samp.num.cumul) ~ treatment + week, data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glmm.nointxn, glmm.notreat, glmm.noweek, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glm.full, test="Chisq") 

# Final Model
Anova(glmm.full, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmm.full)
exp(fixef(glmm.full)) #returns odds ratis instead of logit scale
exp(confint(glmm.full)) #returns odds ratis instead of logit scale
ranef(glmm.full)
VarCorr(glmm.full)



# ANALYZE DATA: Effect of rearing density, larval duration, week, and juv tank id nested within clutch on % surviving individuals weekly ---------------------

# model definition - setting one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16

glmm.full <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment*mean.days.forelimb + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.nointxn <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + mean.days.forelimb + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.notreat <- glmer(live.num/(seed.num-leth.samp.num) ~ mean.days.forelimb + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.nodays <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.noweek <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + mean.days.forelimb + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.null <- glmer(live.num/(seed.num-leth.samp.num) ~ (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glm.full <- glm(live.num/(seed.num-leth.samp.num) ~ treatment*mean.days.forelimb + postmm.week, data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glmm.nointxn, glmm.notreat, glmm.nodays, glmm.noweek, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glm.full, test="Chisq") 

# Final Model
Anova(glmm.nointxn, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmm.nointxn)
exp(fixef(glmm.nointxn)) #returns odds ratis instead of logit scale
exp(confint(glmm.nointxn)) #returns odds ratis instead of logit scale
ranef(glmm.nointxn)
VarCorr(glmm.nointxn)


# ANALYZE DATA: Effect of rearing density, larval duration, week, and juv tank id nested within clutch on % surviving individuals weekly BUT NOW USING CATEGORY FOR LARVAL DURATION ---------------------

# model definition - setting one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16

glmm.full <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + devo.cat*postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.nointxn <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + devo.cat + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.notreat <- glmer(live.num/(seed.num-leth.samp.num) ~ devo.cat + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.nodays <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.noweek <- glmer(live.num/(seed.num-leth.samp.num) ~ treatment + devo.cat + (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glmm.null <- glmer(live.num/(seed.num-leth.samp.num) ~ (1|clutch:juv.tank.id), data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

glm.full <- glm(live.num/(seed.num-leth.samp.num) ~ treatment*devo.cat + postmm.week, data=survi.data.juv[as.numeric(survi.data.juv$postmm.week) > 0,], na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glmm.nointxn, glmm.notreat, glmm.nodays, glmm.noweek, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glm.full, test="Chisq") 

# Final Model
Anova(glmm.nointxn, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmm.nointxn)
exp(fixef(glmm.nointxn)) #returns odds ratis instead of logit scale
exp(confint(glmm.nointxn)) #returns odds ratis instead of logit scale
ranef(glmm.nointxn)
VarCorr(glmm.nointxn)

# create dataframe of predicted values that can be plotted on ggplot later
new.dataframe = expand.grid(prob.survi = seq(from = 0, to = 1, by = 0.1),
                     treatment = c("low density", "high density"),
                     clutch = c("A", "B", "C"),
                     devo.cat = c("early", "mid", "late"),
                     juv.tank.id = unique(survi.data.juv$juv.tank.id),
                     postmm.week =  seq(from = 1, to = 18, by = 1))
new.dataframe$unique.id = paste(new.dataframe$clutch, new.dataframe$juv.tank.id, sep = ":")

# remove random effect ids that don't exist in the model
new.dataframe = new.dataframe[which(new.dataframe$unique.id %in% unique(row.names(ranef(glmm.nointxn)$"clutch:juv.tank.id")) == TRUE),]

#change unique id to be a factor
new.dataframe$unique.id = factor(new.dataframe$unique.id)

#add predictions using model
new.dataframe$predict.mod1 = predict(glmm.nointxn, type = "response", newdata = new.dataframe)

## [-2] drops response from formula
Designmat <- model.matrix(delete.response(terms(glmm.nointxn)),new.dataframe)
predvar <- diag(Designmat %*% vcov(glmm.nointxn) %*% t(Designmat)) 
new.dataframe$SE <- sqrt(predvar) 
new.dataframe$ymin = new.dataframe$predict.mod1 - new.dataframe$SE
new.dataframe$ymin[new.dataframe$ymin < 0] = 0


# PLOT DATASETS: Effect of water level reduction on survival before and at metamorphosis -----------------------
# x-y plot with summarized mean and +/- 1 se for all metrics

# percent survival from tank seeding to tank close-out without clutch included
ggplot(data = survi.data.tad.exp, aes(y=prop.survi*100, x = week, color = water.level.reduc)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) + 
  facet_grid(rows = vars(treatment)) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = water.level.reduc, group = water.level.reduc), show.legend = FALSE) + #show.legend is FALSE so that won't show up on panel plot later
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = water.level.reduc)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=water.level.reduc), show.legend = FALSE) + #show.legend is FALSE so that won't show up on panel plot later
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
  #expand_limits(y = 0) +
  scale_y_continuous(name = "percent survived") +
  scale_x_continuous(name = "pre-metamorphic age (weeks)", breaks = seq(1, 11, by = 1))


# PLOT DATASETS: Effect of rearing density on survival before, at, and after metamorphosis -----------------------
# x-y plot with summarized mean and +/- 1 se for all metrics

# option 1 percent survival from tank seeding to tank close-out without clutch included
plot.survi1 <- ggplot(data = survi.data.tad.exp, aes(y=prop.survi*100, x = week, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment), show.legend = TRUE) + 
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.02),
        legend.justification = c("left", "bottom"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent survived", limits = c(0, 100)) +
  scale_x_continuous(name = "pre-metamorphic age (weeks)", breaks = seq(1, 11, by = 1))

plot.survi2 <- ggplot(data = survi.data.juv[survi.data.juv$postmm.week > 0,], aes(y=prop.survi*100, x = postmm.week, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment), show.legend = TRUE) + 
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.text = element_text(size=10),
        legend.title = element_blank(),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent survived") +
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = seq(1, 18, by = 1))


# PLOT DATASETS: Effect of developmental speed on survival at and after metamorphosis ---------------------

# x-y plot with smoothed binomial curve WITHOUT clutch
plot.survi3 <- ggplot(data = survi.data.juv[survi.data.juv$postmm.week > 0,], aes(y=prop.survi, x = mean.days.forelimb, color = as.factor(postmm.week))) + 
  geom_point(size = 2.5, alpha = 0.7) +
  stat_smooth(method="glm",
              method.args = list(family = binomial), 
              mapping = aes(color=as.factor(postmm.week), 
                            weight = (seed.num-leth.samp.num)), 
              se=F) +
  scale_color_hue(name = "post-metamorphic age", labels = c("1 week", "2 weeks", "3 weeks", "4 weeks", "5 weeks", "6 weeks", "7 weeks", "8 weeks", "9 weeks", "10 weeks", "11 weeks", "12 weeks", "13 weeks", "14 weeks", "15 weeks", "16 weeks", "17 weeks")
  ) +
  theme_bw() +
  theme(legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent survival", limits = c(0,1), breaks = seq(0, 1, by = .25), labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "mean larval duration (days)")

plot.survi3 <- ggplot() + 
  geom_point(size = 2.5, alpha = 0.7, data = survi.data.juv[survi.data.juv$postmm.week > 0,], aes(y=prop.survi, x = postmm.week, color = mean.days.forelimb)) +
  scale_color_gradient(low = "lightgray", high = "black") +
  scale_fill_manual(values = c("lightgray", "darkgray","black")) +
  scale_linetype_manual(values=c(3,2,1)) +
  geom_point(size = 0.1, alpha = 0, data = new.dataframe, aes(y=prob.survi, x = postmm.week)) +
  geom_smooth(color = "black",
              data = new.dataframe,
              aes(x = postmm.week, y = predict.mod1, linetype=devo.cat),
              se = F) +
  geom_ribbon(data = new.dataframe, alpha = 0.1, aes(x = postmm.week, ymin=ymin, ymax=predict.mod1+SE, fill = devo.cat)) +
  theme_bw() +
  labs(linetype="larval duration category", color="larval duration (days)") +
  theme(legend.position = c(0.98, 0.96),
        legend.justification = c("right", "top"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10),
        legend.spacing.y = unit(0.05, "cm"), 
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent survived", limits = c(0,1), breaks = seq(0, 1, by = .25), labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = seq(1, 18, by = 1))



# PLOT DATA: Create panel plot with all survival data across all sampling points --------------
ggarrange(plot.survi1, plot.survi2, plot.survi3,
          ncol = 3,
          nrow = 1,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b", "c"),
          font.label = list(size = 20, color = "black"))
