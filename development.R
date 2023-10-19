# this script creates figures of the developmental results for the 2023 ranid developmental paper
# includes effects of rearing density on developmental timepoints (1. time to forelimb emergence, time to complete metamorphosis, laterality of forelimb emergence)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for developmental data and survival data (provides tank-level percent metamorphosed by end of experiment)
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.exp = read.csv("Database_Survivorship - Tadpole Experiment Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.tad = read.csv("Database_Survivorship - Tadpole Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
devo.data = devo.data[devo.data$gs.code == "RS",]
survi.data.exp = survi.data.exp[survi.data.exp$gs.code == "RS",]
survi.data.tad = survi.data.tad[survi.data.tad$gs.code == "RS",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data$treatment[devo.data$treatment == "control"] = "low density"
survi.data.exp$treatment[survi.data.exp$treatment == "control"] = "low density"
survi.data.tad$treatment[survi.data.tad$treatment == "control"] = "low density"

#change column classes
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)
devo.data$days.forelimb = as.integer(devo.data$days.forelimb)
devo.data$days.forelimb.tail = as.integer(devo.data$days.forelimb.tail)

survi.data.exp$treatment = factor(survi.data.exp$treatment)
survi.data.exp$larv.tank.id = factor(survi.data.exp$larv.tank.id)
survi.data.exp$prop.surv.forelimb = as.numeric(survi.data.exp$prop.surv.forelimb)

survi.data.tad$treatment = factor(survi.data.tad$treatment)
survi.data.tad$larv.tank.id = factor(survi.data.tad$larv.tank.id)
survi.data.tad$prop.surv.forelimb = as.numeric(survi.data.tad$prop.surv.forelimb)

# create unique tank id for survi.data.exp and devo.data
survi.data.exp$unique.id = paste(survi.data.exp$gs.code, survi.data.exp$clutch, survi.data.exp$larv.tank.id, survi.data.exp$week, sep = "_")
devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
survi.data.tad$unique.id = paste(survi.data.tad$gs.code, survi.data.tad$clutch, survi.data.tad$larv.tank.id, survi.data.tad$week, sep = "_")

# rename two columns to match survi.data.tad for later binding
colnames(survi.data.exp)[11] = "photo.num"
colnames(survi.data.exp)[12] = "leth.samp.num.cumul"



# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to survi.data.tad.exp -----------------------

# create column to store the developmental data, mean.days.forelimb represents days.forelimb for each individual (rather than an average, since we actually have individual level data for this) but need to keep column name the same as in the juvenile dataset so we can combine them in the next step
survi.data.exp$mean.days.forelimb = NA

# fill developmental data column for survi.data.tad.exp
for(i in 1:length(unique(survi.data.exp$unique.id))){
  survi.data.exp$mean.days.forelimb[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(survi.data.exp$unique.id)[i] & devo.data$first.six == "yes"], na.rm = TRUE)
}


# COMPILE DATA: Generate end of experiment metamorphosis dataframe in terms of 0 (dead) and 1 (alive) -----------------------
#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.exp[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[1],]
temp3 = data.frame(gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                   gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                   clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                   treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                   treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                   larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num),
                   mean.days.forelimb = rep(temp$mean.days.forelimb, temp$seed.num - temp$leth.samp.num),
                   status = NA
)
temp3$status[0:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.exp$unique.id))){
  temp = survi.data.exp[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i],]
  temp2 = data.frame(gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num),
                     mean.days.forelimb = rep(temp$mean.days.forelimb, temp$seed.num - temp$leth.samp.num),
                     status = NA
  )
  temp2$status[0:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.exp.longform = temp3
rm(temp3)


# COMPILE DATA: Generate end of experiment metamorphosis dataframe in terms of 0 (dead) and 1 (alive) BUT ONLY CONSIDERING SURVIVING INDIVIDUALS OUT OF SEEDED -----------------------
#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.exp[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[1],]
temp3 = data.frame(gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                   gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                   clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                   treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                   mean.days.forelimb = rep(temp$mean.days.forelimb, temp$photo.num + temp$metamorph.num.cumul),
                   status = NA
)
temp3$status[1:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.exp$unique.id))){
  temp = survi.data.exp[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i],]
  temp2 = data.frame(gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                     gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                     clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                     treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                     mean.days.forelimb = rep(temp$mean.days.forelimb, temp$photo.num + temp$metamorph.num.cumul),
                     status = NA
  )
  temp2$status[1:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.exp.longformsurvi = temp3
rm(temp3)


# COMPILE DATA: Generate weekly metamorphosis dataframe in terms of 0 (dead) and 1 (alive) -----------------------

#combine weekly data with end of experiment data
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
                                               "unique.id")])

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (metamorphosed) is calculated from the metamorphosed individuals
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
temp3$status[0:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                    gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num),
                     status = NA
  )
  temp2$status[0:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.tad.exp.longform = temp3
rm(temp3)


# COMPILE DATA: Generate weekly metamorphosis dataframe in terms of 0 (dead) and 1 (alive) BUT ONLY CONSIDERING SURVIVING INDIVIDUALS OUT OF SEEDED -----------------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[1],]
temp3 = data.frame(week = rep(temp$week, temp$photo.num + temp$metamorph.num.cumul),
                   gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                   gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                   clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                   treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                   status = NA
)
temp3$status[0:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.tad.exp$unique.id)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$photo.num + temp$metamorph.num.cumul),
                     gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                     gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                     clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                     treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                     status = NA
  )
  temp2$status[0:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.tad.exp.longformsurvi = temp3
rm(temp3)



# ANALYZE DATA: Effect of rearing density and clutch:larval tank on larval duration ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but ML does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
lmer.full <- lmer(days.forelimb ~ treatment + (1|clutch:larv.tank.id), data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)
lmer.null<- lmer(days.forelimb ~ (1|clutch:larv.tank.id), data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)
lm.null <- lm(days.forelimb ~ treatment, data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmer.full, lmer.null, test="Chisq")

# checking whether random effect important to include
anova(lmer.full, lm.null, test = "Chisq")

# Check Model Assumptions
check_model(lmer.full)

# Final Model
summary(lmer.full)
Anova(lmer.full, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
fixef(lmer.full)
ranef(lmer.full)
confint(lmer.full)
get_variance(lmer.full)


# ANALYZE DATA: Effect of rearing density and clutch:larval tank on metamorphosis duration ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but ML does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
lmer.full <- lmer(days.forelimb.tail ~ treatment + (1|clutch:larv.tank.id), data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)
lmer.null<- lmer(days.forelimb.tail ~ (1|clutch:larv.tank.id), data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)
lm.null <- lm(days.forelimb.tail ~ treatment, data = devo.data[devo.data$first.six == "yes",], na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test + checking whether random effect important to include
anova(lmer.full, lmer.null, test="Chisq")
anova(lmer.full, lm.null, test = "Chisq")

# Check Model Assumptions
check_model(lmer.full)

# Final Model
summary(lmer.full)
Anova(lmer.full, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
fixef(lmer.full)
ranef(lmer.full)
confint(lmer.full)
get_variance(lmer.full)

summary(lmer.null)
Anova(lmer.null, type = "II") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
fixef(lmer.null)
ranef(lmer.null)
confint(lmer.null)
get_variance(lmer.null)


# ANALYZE DATA: Effect of rearing density, larval duration, week, and larval tank id nested within clutch on likelihood of metamorphosing with forelimb emerged weekly and by end of experiment ---------------------

# model definition - setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
cox.model.mixed.full <- coxme(Surv(status) ~ treatment + week + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model.mixed.noweek <- coxme(Surv(status) ~ treatment + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model.mixed.notreat <- coxme(Surv(status) ~ week + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model.mixed.null <- coxme(Surv(status) ~ (1|clutch/larv.tank.id), data=survi.data.tad.exp.longform)
cox.model <- coxph(Surv(status) ~ treatment + week, data=survi.data.tad.exp.longform)

# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model.mixed.noweek,cox.model.mixed.notreat, cox.model.mixed.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model, test="Chisq") 

# Final Model
Anova(cox.model.mixed.full, type = "II")
summary(cox.model.mixed.full)
exp(confint(cox.model.mixed.full))
VarCorr(cox.model.mixed.full)


# ANALYZE DATA: Effect of rearing density, larval duration, week, and larval tank id nested within clutch on likelihood of metamorphosing based on seeded AND SURVIVED individuals with forelimb emerged weekly and by end of experiment ---------------------

# model definition - setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.
cox.model.mixed.full <- coxme(Surv(status) ~ treatment + week + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longformsurvi)
cox.model.mixed.noweek <- coxme(Surv(status) ~ treatment + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longformsurvi)
cox.model.mixed.notreat <- coxme(Surv(status) ~ week + (1|clutch/larv.tank.id), data=survi.data.tad.exp.longformsurvi)
cox.model.mixed.null <- coxme(Surv(status) ~ (1|clutch/larv.tank.id), data=survi.data.tad.exp.longformsurvi)
cox.model <- coxph(Surv(status) ~ treatment + week, data=survi.data.tad.exp.longformsurvi)

# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model.mixed.notank, cox.model.mixed.noweek, cox.model.mixed.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(cox.model.mixed.full, cox.model, test="Chisq") 

# Final Model
Anova(cox.model.mixed.full, type = "II")
summary(cox.model.mixed.full)
exp(confint(cox.model.mixed.full))
VarCorr(cox.model.mixed.full)


# ANALYZE DATA: Effect of rearing density, larval duration, week, and larv tank id nested within clutch on % metamorphosing individuals weekly up until metamorphosis ---------------------

# model definition - setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16
glmm.full <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment + week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.noweek <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.notreat <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.null <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glm.full <- glm(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week, data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glmm.nointxn, glmm.noweek, glmm.notreat, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glm.full, test="Chisq") 

# Final Model
Anova(glmm.full, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
summary(glmm.full)
exp(fixef(glmm.full)) #returns odds ratios instead of logit scale
exp(confint(glmm.full)) #returns odds ratios instead of logit scale
ranef(glmm.full)
VarCorr(glmm.full)



# ANALYZE DATA: Effect of rearing density, larval duration, week, and larv tank id nested within clutch on % metamorphosing individuals CALCULATED ONLY FROM SURVIVING INDIVIDUALS weekly up until metamorphosis ---------------------

# model definition - setting one random effect as week and one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16
glmm.full <- glmer(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ treatment*week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))

glmm.nointxn <- glmer(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ treatment + week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))

glmm.noweek <- glmer(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ treatment + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))

glmm.notreat <- glmer(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ week + (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))

glmm.null <- glmer(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ (1|clutch:larv.tank.id), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))

glm.full <- glm(metamorph.num.cumul/(photo.num+metamorph.num.cumul) ~ treatment*week, data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (photo.num+metamorph.num.cumul))


# model selection using likelihood ratio test; higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glmm.nointxn, glmm.noweek, glmm.notreat, glmm.null, test="Chisq")

# check if random effect significant - if model without the random effect is not significantly different than model with the random effect, then random effect isn't explaining the variance much. Even if not significant, should include in final model to more accurately account for covariance. higher negative log-lik (closer to 0) provides better fit to the data 
anova(glmm.full, glm.full, test="Chisq") 

# Final Model
Anova(glmm.full, type = "II")
summary(glmm.full)
exp(fixef(glmm.full)) #returns odds ratios instead of logit scale
exp(confint(glmm.full)) #returns odds ratios instead of logit scale
ranef(glmm.full)
VarCorr(glmm.full)




# PLOT DATASETS: Effect of rearing density on larval duration -----------------------

# plotted with clutch separated under each density group 
plot.devo.1 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = clutch, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days)") +
  scale_x_discrete(name = "clutches separated")

# plotted with clutch lumped together under each density group 
plot.devo.2 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = treatment, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=16, color = "white"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days)") +
  scale_x_discrete(name = "clutches combined")


# PLOT DATASETS: Effect of rearing density on proportion tank initated metamorphosis by experiment completion ---------------------
# x-y plot with summarized mean and +/- 1 se for all metrics

# percent metamorphosed from tank seeding to tank close-out without clutch included
plot.devo.3 <- ggplot(data = survi.data.tad.exp, aes(y=prop.seed.forelimb*100, x = week, color = treatment)) + 
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
        axis.text.x=element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "percent metamorphosed", limits = c(0,30)) +
  scale_x_continuous(name = "age (weeks)", breaks = seq(1, 11, by = 1))

# raw number metamorphosed from tank seeding to tank close-out without clutch included
plot.devo.4 <- ggplot(data = survi.data.tad.exp, aes(y=metamorph.num.cumul, x = week, color = treatment)) + 
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
        axis.text.x=element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "number metamorphosed", limits = c(0,30)) +
  scale_x_continuous(name = "age (weeks)", breaks = seq(1, 11, by = 1))


# PLOT DATASETS: Panel plot of little-to-no effect of rearing density on developmental metrics ---------------------
ggarrange(plot.devo.3, plot.devo.4, plot.devo.1, plot.devo.2,
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))
# show.legend = FALSE in boxplots so that the common legend pulls from the other plots



# SUMMARY TABLES: Developmental metrics ---------------------

#reset working directory to be outputs folder
setwd("~/Desktop/R Working Directory/Outputs")

devo.data.summ = devo.data %>%
  group_by(juv.tank.id) %>%
  filter(first.six == "yes" & is.na(juv.tank.id)==FALSE) %>%
  summarise(mean.days.forelimb = round(mean(days.forelimb, na.rm = TRUE), 2),
            sd.days.forelimb = round(sd(days.forelimb, na.rm = TRUE), 2),
            min.days.forelimb = round(min(days.forelimb, na.rm = TRUE), 2),
            max.days.forelimb = round(max(days.forelimb, na.rm = TRUE), 2),
            range.days.forelimb = round(max(days.forelimb, na.rm = TRUE), 2) - round(min(days.forelimb, na.rm = TRUE), 2),
            mean.days.completion = round(mean(days.forelimb.tail, na.rm = TRUE), 2),
            sd.days.completion = round(sd(days.forelimb.tail, na.rm = TRUE), 2)
  )

survi.data.exp.summ = survi.data.exp %>%
  group_by(treatment) %>%
  summarise(num.seeded = sum(seed.num - leth.samp.num, na.rm = TRUE), 
            num.survi = sum(larv.num + metamorph.num.cumul, na.rm = TRUE),
            num.metamorph = sum(metamorph.num.cumul, na.rm = TRUE),
            mean.metamorph = round(mean(metamorph.num.cumul, na.rm = TRUE), 2),
            sd.metamorph = round(sd(metamorph.num.cumul, na.rm = TRUE), 2),
            mean.prop.seeded = round(mean(prop.seed.forelimb, na.rm = TRUE), 3),
            sd.prop.seeded = round(sd(prop.seed.forelimb, na.rm = TRUE), 3),
            mean.prop.survi = round(mean(prop.surv.forelimb, na.rm = TRUE), 3),
            sd.prop.survi = round(sd(prop.surv.forelimb, na.rm = TRUE), 3)
  )