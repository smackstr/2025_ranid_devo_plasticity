rm(list = ls())

# LOAD PACKAGES & SET WD ---------------------

requiredPackages = c("lme4", "car", "DHARMa", "MuMIn", "tidyverse", "purrr")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}

#set working directory to wherever databases are stored
wd = #add filepath here in ""
setwd(wd)

# LOAD DATABASES ---------------------

# read in cleaned/joined databases
morph.data.tad = read.csv("Database_Morphometrics_Tadpole.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.mm = read.csv("Database_Morphometrics_Metamorph.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv = read.csv("Database_Morphometrics_Juvenile.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.mm.juv = read.csv("Database_Morphometrics_MetamorphAndJuvenile.csv", header = TRUE, skip = 0, na.strings = "NA")

devo.data = read.csv("Database_Development.csv", header = TRUE, skip = 0, na.strings = "NA")

survi.data.tad = read.csv("Database_Survivorship_Tadpole_Weekly.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.exp = read.csv("Database_Survivorship_Tadpole_Experiment.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.tad.exp = read.csv("Database_Survivorship_Tadpole_WeeklyAndExperiment.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.tad.exp.longform = read.csv("Database_Survivorship_Tadpole_WeeklyAndExperiment_Longform.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.tad.exp.longformsurvi = read.csv("Database_Survivorship_Tadpole_WeeklyAndExperiment_LongformSurvi.csv", header = TRUE, skip = 0, na.strings = "NA")

survi.data.juv = read.csv("Database_Survivorship_Juvenile_Weekly.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.juv.longform = read.csv("Database_Survivorship_Juvenile_Longform.csv", header = TRUE, skip = 0, na.strings = "NA")

# ADJUST COLUMN CLASSES so reference level is consistent across analyses ---------------------

morph.data.tad <- morph.data.tad %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

morph.data.mm <- morph.data.mm %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

morph.data.juv <- morph.data.juv %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

morph.data.mm.juv <- morph.data.mm.juv %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.tad <- survi.data.tad %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.exp <- survi.data.exp %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.tad.exp <- survi.data.tad.exp %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.tad.exp.longform <- survi.data.tad.exp.longform %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.tad.exp.longformsurvi <- survi.data.tad.exp.longformsurvi %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.juv <- survi.data.juv %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

survi.data.juv.longform <- survi.data.juv.longform %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density"))

devo.data <- devo.data %>%
  mutate_at(vars(clutch, clutchtank), factor) %>%
  mutate_at(vars(treatment), factor, levels = c("low density", "high density")) %>%
  mutate_at(vars(first.six, first.10perc), factor, levels = c("yes", "no"))

# Supp Table 1a. ANALYZE DATA - LARVAL DEVELOPMENT : Effect of rearing density on metamorphosis status (0 = tadpole, 1 = metamorphosis initiated) of seeded individuals of each larval tank ---------------------

#define candidate models
glmm.full.01 <- glmer(status ~ treatment*week + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes <- glmer(status ~ treatment*week + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + week  + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

glmm.null.01 <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc 
model.sel = arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
             glmm.nointxn.01, glmm.nointxn.01.slopes,
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

#check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tad.exp.longform$treatment)
plotResiduals(glmm.full.01, form = survi.data.tad.exp.longform$week)
testDispersion(glmm.full.01)

# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant


# Supp Table 1b. ANALYZE DATA - LARVAL DEVELOPMENT : Effect of rearing density on metamorphosis status (0 = tadpole, 1 = metamorphosis initiated) of surviving individuals of each larval tank ---------------------

#define candidate models
glmm.full.01 <- glmer(status ~ treatment*week + (1|clutchtank), data=survi.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes <- glmer(status ~ treatment*week + (week||clutchtank), data=survi.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + week + (1|clutchtank), data=survi.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + (week||clutchtank), data=survi.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.01 <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc 
model.sel = arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
             glmm.nointxn.01, glmm.nointxn.01.slopes,
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

#check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(final.mod, form = survi.data.tad.exp.longformsurvi$treatment[is.na(survi.data.tad.exp.longformsurvi$treatment)==FALSE])
testDispersion(final.mod)

# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant


# Supp Table 2. ANALYZE DATA - LARVAL DEVELOPMENT : Effect of rearing density on larval duration for first 10% of each larval tank to metamorphose ---------------------

#define candidate models
lmm.full <- lmer(days.forelimb ~ treatment + (1|clutchtank), data = filter(devo.data, first.10perc == "yes"), na.action = na.omit)
lmm.null <- lmer(days.forelimb ~ 1 + (1|clutchtank), data = filter(devo.data, first.10perc == "yes"), na.action = na.omit)

# nested model selection using likelihood-ratio test
model.sel = arrange(anova(lmm.full, lmm.null), AIC)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(final.mod) #tests for over- and under-dispersion
testZeroInflation(final.mod) #tests if more zeroes than expected
testCategorical(final.mod, catPred = devo.data$treatment)
testOutliers(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 3. ANALYZE DATA - LARVAL DEVELOPMENT: Effect of rearing density and larval duration on tank-level average developmental stage (Gosner) at end of experiment ---------------------

#define candidate models
lmm.full <- lmer(tank.devo.stage ~ treatment + (1|clutch), data=survi.data.exp, na.action = na.omit)

lmm.null <- lmer(tank.devo.stage ~ 1 + (1|clutch), data=survi.data.exp, na.action = na.omit)

# nested model selection using likelihood-ratio test
model.sel = arrange(anova(lmm.full, lmm.null), AIC)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(final.mod) #tests for over- and under-dispersion
testZeroInflation(final.mod) #tests if more zeroes than expected
testCategorical(final.mod, catPred = survi.data.exp$treatment)
testOutliers(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 4a. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on TOTAL LENGTH ---------------------

#define candidate models
lmm.full <- lmer(tl.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(tl.cm ~ treatment*week + (week||clutch), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #singular fit with more complex random effect structure, so simplified

lmm.full.log <- lmer(log(tl.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(tl.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.log <- lmer(log(tl.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(tl.cm) ~ treatment*log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lmm.full.log.poly <- lmer(log(tl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(tl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx <- lmer(tl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(tl.cm ~ treatment + week + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex random effect structure, so simplified

lmm.nointx.log <- lmer(log(tl.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(tl.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.log <- lmer(log(tl.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(tl.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.poly <- lmer(log(tl.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(tl.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.null <- lmer(tl.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(tl.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes,
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes,
             lmm.full.log.poly, lmm.full.log.poly.slopes,
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.slopes, quantreg=T, plot = T)
testDispersion(lmm.nointx.slopes)
testZeroInflation(lmm.nointx.slopes)
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$treatment) 
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$week)

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly)
testZeroInflation(lmm.full.log.poly)
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$treatment) 
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$week)

final.mod = lmm.full.log.poly

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 4b. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on BODY LENGTH ---------------------

#define candidate models
lmm.full <- lmer(bl.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(bl.cm ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log <- lmer(log(bl.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(bl.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.log <- lmer(log(bl.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(bl.cm) ~ treatment*log(week) + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex random effect structure, so simplified

lmm.full.log.poly <- lmer(log(bl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(bl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex AND simplified random effect structure

lmm.nointx <- lmer(bl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(bl.cm ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log <- lmer(log(bl.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(bl.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.log <- lmer(log(bl.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(bl.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.poly <- lmer(log(bl.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(bl.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(bl.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(bl.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes, lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes,
             lmm.full.log.poly, 
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.slopes, quantreg=T, plot = T)
testDispersion(lmm.nointx.slopes)
testZeroInflation(lmm.nointx.slopes)
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$treatment) 
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$week)

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly)
testZeroInflation(lmm.full.log.poly)
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$treatment) 
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$week)

final.mod = lmm.full.log.poly

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 4c. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on TAIL LENGTH ---------------------

#define candidate models
lmm.full <- lmer(tail.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(tail.cm ~ treatment*week + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex random effect structure, so simplified

lmm.full.log <- lmer(log(tail.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(tail.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.full.log.log <- lmer(log(tail.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(tail.cm) ~ treatment*log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.poly <- lmer(log(tail.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(tail.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.nointx <- lmer(tail.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(tail.cm ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log <- lmer(log(tail.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(tail.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.log <- lmer(log(tail.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(tail.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.nointx.log.poly <- lmer(log(tail.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(tail.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(tail.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(tail.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full,lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes,
             lmm.full.log.poly, lmm.full.log.poly.slopes,
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.slopes, quantreg=T, plot = T)
testDispersion(lmm.nointx.slopes)
testZeroInflation(lmm.nointx.slopes)
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$treatment) 
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$week)

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly.slopes, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly.slopes)
testZeroInflation(lmm.full.log.poly.slopes)
testCategorical(lmm.full.log.poly.slopes, catPred = morph.data.tad$treatment) 
testCategorical(lmm.full.log.poly.slopes, catPred = morph.data.tad$week)

final.mod = lmm.full.log.poly.slopes

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 4d. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on HEAD WIDTH ---------------------

#define candidate models
lmm.full <- lmer(hw.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(hw.cm ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log <- lmer(log(hw.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(hw.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.log <- lmer(log(hw.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(hw.cm) ~ treatment*log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.poly <- lmer(log(hw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(hw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx <- lmer(hw.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(hw.cm ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log <- lmer(log(hw.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(hw.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.log <- lmer(log(hw.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(hw.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.poly <- lmer(log(hw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(hw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(hw.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(hw.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc
arrange(AICc(lmm.full, lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes,
             lmm.full.log.poly, lmm.full.log.poly.slopes,
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.slopes, quantreg=T, plot = T)
testDispersion(lmm.nointx.slopes)
testZeroInflation(lmm.nointx.slopes)
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$treatment) 
testCategorical(lmm.nointx.slopes, catPred = morph.data.tad$week)

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly)
testZeroInflation(lmm.full.log.poly)
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$treatment) 
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$week)

final.mod = lmm.full.log.poly

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 4e. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on BODY WIDTH ---------------------

#define candidate models
lmm.full <- lmer(bw.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(bw.cm ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.full.log <- lmer(log(bw.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(bw.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.full.log.log <- lmer(log(bw.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(bw.cm) ~ treatment*log(week) + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex random effect structure, so simplified

lmm.full.log.poly <- lmer(log(bw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(bw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutch), data = morph.data.tad, na.action = na.omit) #singular fit with more complex random effect structure, so simplified

lmm.nointx <- lmer(bw.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(bw.cm ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log <- lmer(log(bw.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(bw.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.nointx.log.log <- lmer(log(bw.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(bw.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.poly <- lmer(log(bw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(bw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(bw.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(bw.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc non-transformed response variable
arrange(AICc(lmm.full, lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes, 
             lmm.full.log.poly, lmm.full.log.poly.slopes,
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T)
testDispersion(lmm.nointx)
testZeroInflation(lmm.nointx)
testCategorical(lmm.nointx, catPred = morph.data.tad$treatment) 
testCategorical(lmm.nointx, catPred = morph.data.tad$week)

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly)
testZeroInflation(lmm.full.log.poly)
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$treatment) 
testCategorical(lmm.full.log.poly, catPred = morph.data.tad$week)

final.mod = lmm.full.log.poly

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 4f. ANALYZE DATA - TADPOLE MORPHOMETRICS: Effect of rearing density and age on TAIL-MUSCLE WIDTH ---------------------

#define candidate models
lmm.full <- lmer(tmw.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.slopes <- lmer(tmw.cm ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log <- lmer(log(tmw.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.slopes <- lmer(log(tmw.cm) ~ treatment*week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.log <- lmer(log(tmw.cm) ~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log.slopes <- lmer(log(tmw.cm) ~ treatment*log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.full.log.poly <- lmer(log(tmw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(tmw.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx <- lmer(tmw.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.slopes <- lmer(tmw.cm ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log <- lmer(log(tmw.cm) ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.slopes <- lmer(log(tmw.cm) ~ treatment + week + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.log <- lmer(log(tmw.cm) ~ treatment + log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.log.slopes <- lmer(log(tmw.cm) ~ treatment + log(week) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx.log.poly <- lmer(log(tmw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.nointx.log.poly.slopes <- lmer(log(tmw.cm)~ treatment + poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit) 

lmm.null <- lmer(tmw.cm ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.null.log <- lmer(log(tmw.cm) ~ 1 + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.slopes,
             lmm.nointx, lmm.nointx.slopes, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.slopes,
             lmm.full.log.log, lmm.full.log.log.slopes,
             lmm.full.log.poly, lmm.full.log.poly.slopes,
             lmm.nointx.log, lmm.nointx.log.slopes,
             lmm.nointx.log.log, lmm.nointx.log.log.slopes,
             lmm.nointx.log.poly, lmm.nointx.log.poly.slopes,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.slopes, quantreg=T, plot = T)
testDispersion(lmm.nointx.slopes)
testZeroInflation(lmm.nointx.slopes)
testCategorical(lmm.nointx.slopes, catPred =  morph.data.tad$treatment[is.na(morph.data.tad$tmw.cm) == FALSE]) 
testCategorical(lmm.nointx.slopes, catPred =  morph.data.tad$week[is.na(morph.data.tad$tmw.cm) == FALSE]) 

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.full.log.poly.slopes, quantreg=T, plot = T)
testDispersion(lmm.full.log.poly.slopes)
testZeroInflation(lmm.full.log.poly.slopes)
testCategorical(lmm.full.log.poly.slopes, catPred = morph.data.tad$treatment[is.na(morph.data.tad$tmw.cm) == FALSE]) 
testCategorical(lmm.full.log.poly.slopes, catPred =  morph.data.tad$week[is.na(morph.data.tad$tmw.cm) == FALSE]) 

final.mod = lmm.full.log.poly.slopes

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 5. ANALYZE DATA - LARVAL DEVELOPMENT : Effect of rearing density on larval duration for first six individuals of each larval tank to metamorphose ---------------------

#define candidate models
lmm.full <- lmer(days.forelimb ~ treatment + (1|clutchtank), data = filter(devo.data, first.six == "yes"), na.action = na.omit)
lmm.null <- lmer(days.forelimb ~ 1 + (1|clutchtank), data = filter(devo.data, first.six == "yes"), na.action = na.omit)

# nested model selection using likelihood-ratio test
anova(lmm.full, lmm.null)

final.mod = lmm.full #parameterized model not significantly better than null model, but listing final.mod to report data in Supp. Table 5

# check assumptions
simulateResiduals(final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(final.mod) #tests for over- and under-dispersion
testZeroInflation(final.mod) #tests if more zeroes than expected
testOutliers(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 6a. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density and larval duration on MASS ---------------------

#define candidate models
lmm.full <- lmer(mass.g ~ treatment*days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(mass.g ~ treatment*days.forelimb + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.full.log <- lmer(log(mass.g) ~ treatment*days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.log.slopes <- lmer(log(mass.g) ~ treatment*days.forelimb + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointx <- lmer(mass.g ~ treatment + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointx.slopes <- lmer(mass.g ~ treatment + days.forelimb + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointx.log <- lmer(log(mass.g) ~ treatment + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointx.log.slopes <- lmer(log(mass.g) ~ treatment + days.forelimb + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.null<- lmer(mass.g ~ (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.null.log <- lmer(log(mass.g) ~ (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, 
             lmm.nointx, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, 
             lmm.nointx.log,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T)
testDispersion(lmm.nointx)
testZeroInflation(lmm.nointx)
testCategorical(lmm.nointx, catPred = morph.data.mm$treatment[is.na(morph.data.mm$mass.g) == FALSE]) 

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.log, quantreg=T, plot = T)
testDispersion(lmm.nointx.log)
testZeroInflation(lmm.nointx.log)
testCategorical(lmm.nointx.log, catPred = morph.data.mm$treatment[is.na(morph.data.mm$mass.g) == FALSE]) 

final.mod = lmm.nointx

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 6b. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density, larval duration, and metamorphic mass on BODY CONDITION (SCALED MASS INDEX, SMI) ---------------------

#define candidate models
lmm.full <- lmer(smi ~ treatment*days.forelimb*mass.g +  (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(smi ~ treatment*days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn1 <- lmer(smi ~ treatment*days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn1.slopes <- lmer(smi ~ treatment*days.forelimb + mass.g + (treatment|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn2 <- lmer(smi ~ treatment + days.forelimb*mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn2.slopes <- lmer(smi ~ treatment + days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn3 <- lmer(smi ~ treatment*mass.g + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn3.slopes <- lmer(smi ~ treatment*mass.g + days.forelimb + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn <- lmer(smi ~ treatment + days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.slopes <- lmer(smi ~ treatment + days.forelimb + mass.g + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.null <- lmer(smi ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc for non-transformed response variable
model.sel = arrange(AICc(lmm.full, 
             lmm.nointxn, lmm.nointxn1, lmm.nointxn2, lmm.nointxn3,
             lmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = morph.data.mm$treatment[is.na(morph.data.mm$smi) == FALSE]) 

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 6c. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density, larval duration, and metamorphic mass on SNOUT-VENT LENGTH ---------------------

#define candidate models
lmm.full <- lmer(svl.mm ~ treatment*days.forelimb*mass.g +  (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(svl.mm ~ treatment*days.forelimb*mass.g +  (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn1 <- lmer(svl.mm ~ treatment*days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn1.slopes <- lmer(svl.mm ~ treatment*days.forelimb + mass.g + (treatment|clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn2 <- lmer(svl.mm ~ treatment + days.forelimb*mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn2.slopes <- lmer(svl.mm ~ treatment + days.forelimb*mass.g +  (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn3 <- lmer(svl.mm ~ treatment*mass.g + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn3.slopes <- lmer(svl.mm ~ treatment*mass.g + days.forelimb + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.nointxn <- lmer(svl.mm ~ treatment + days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.slopes <- lmer(svl.mm ~ treatment + days.forelimb + mass.g + (treatment||clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE) #singular fit with more complex AND simplified random effect

lmm.null <- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc for non-transformed response variable
model.sel = arrange(AICc(lmm.full, 
                         lmm.nointxn, lmm.nointxn1, lmm.nointxn2, lmm.nointxn3,
                         lmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)
testCategorical(final.mod, catPred = morph.data.mm$treatment[is.na(morph.data.mm$svl.mm) == FALSE]) 

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 6d. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density, larval duration, and metamorphic mass on FOREARM LENGTH ---------------------

#define candidate models
lmm.full <- lmer(r.forelimb.mm ~ treatment*days.forelimb*mass.g +  (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(r.forelimb.mm ~ treatment*days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #failure to converge

lmm.nointxn1 <- lmer(r.forelimb.mm ~ treatment*days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn1.slopes <- lmer(r.forelimb.mm ~ treatment*days.forelimb + mass.g + (treatment|clutch), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #singular fit with more complex and simplified random effect structure

lmm.nointxn2 <- lmer(r.forelimb.mm ~ treatment + days.forelimb*mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn2.slopes <- lmer(r.forelimb.mm ~ treatment + days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lmm.nointxn3 <- lmer(r.forelimb.mm ~ treatment*mass.g + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn3.slopes <- lmer(r.forelimb.mm ~ treatment*mass.g + days.forelimb + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lmm.nointxn <- lmer(r.forelimb.mm ~ treatment + days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.slopes <- lmer(r.forelimb.mm ~ treatment + days.forelimb + mass.g + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lmm.null <- lmer(r.forelimb.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc 
model.sel = arrange(AICc(lmm.full, lmm.full.slopes,
                         lmm.nointxn1,
                         lmm.nointxn2, lmm.nointxn2.slopes,
                         lmm.nointxn3, lmm.nointxn3.slopes,
                         lmm.nointxn, lmm.nointxn.slopes,
                         lmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)


# Supp Table 6e. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density, larval duration, and metamorphic mass on RIGHT TIBIA LENGTH ---------------------

#define candidate models
lmm.full <- lmer(r.tibia.mm ~ treatment*days.forelimb*mass.g +  (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(r.tibia.mm ~ treatment*days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #failure to converge

lmm.nointxn1 <- lmer(r.tibia.mm ~ treatment*days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn1.slopes <- lmer(r.tibia.mm ~ treatment*days.forelimb + mass.g + (treatment|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn2 <- lmer(r.tibia.mm ~ treatment + days.forelimb*mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn2.slopes <- lmer(r.tibia.mm ~ treatment + days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.nointxn3 <- lmer(r.tibia.mm ~ treatment*mass.g + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn3.slopes <- lmer(r.tibia.mm ~ treatment*mass.g + days.forelimb + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(r.tibia.mm ~ treatment + days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.slopes <- lmer(r.tibia.mm ~ treatment + days.forelimb + mass.g + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.null <- lmer(r.tibia.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc 
model.sel = arrange(AICc(lmm.full,
                         lmm.nointxn1, lmm.nointxn1.slopes,
                         lmm.nointxn2, lmm.nointxn2.slopes,
                         lmm.nointxn3, lmm.nointxn3.slopes,
                         lmm.nointxn, lmm.nointxn.slopes,
                         lmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "III")  #type 3 because interaction term significant
summary(final.mod)


# Supp Table 6f. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of rearing density, larval duration, and metamorphic mass on RIGHT THIGH LENGTH ---------------------

#define candidate models
lmm.full <- lmer(r.thigh.mm ~ treatment*days.forelimb*mass.g +  (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.slopes <- lmer(r.thigh.mm ~ treatment*days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

lmm.nointxn1 <- lmer(r.thigh.mm ~ treatment*days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn1.slopes <- lmer(r.thigh.mm ~ treatment*days.forelimb + mass.g + (treatment|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn2 <- lmer(r.thigh.mm ~ treatment + days.forelimb*mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn2.slopes <- lmer(r.thigh.mm ~ treatment + days.forelimb*mass.g +  (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn3 <- lmer(r.thigh.mm ~ treatment*mass.g + days.forelimb + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn3.slopes <- lmer(r.thigh.mm ~ treatment*mass.g + days.forelimb + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(r.thigh.mm ~ treatment + days.forelimb + mass.g + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.slopes <- lmer(r.thigh.mm ~ treatment + days.forelimb + mass.g + (treatment||clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 

lmm.null <- lmer(r.thigh.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc 
model.sel = arrange(AICc(lmm.full, lmm.full.slopes,
                         lmm.nointxn1, lmm.nointxn1.slopes,
                         lmm.nointxn2, lmm.nointxn2.slopes,
                         lmm.nointxn3, lmm.nointxn3.slopes,
                         lmm.nointxn, lmm.nointxn.slopes,
                         lmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T)
testDispersion(final.mod)
testZeroInflation(final.mod)

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 7. ANALYZE DATA - METAMORPH MORPHOMETRICS: Effect of water level reduction and larval duration on MASS for low-density tanks (since high density did not experience water level reduction)---------------------

#define candidate models
lmm.full <- lmer(mass.g ~ days.forelimb*water.level.reduc + (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)
lmm.full.log <- lmer(log(mass.g) ~ days.forelimb*water.level.reduc + (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)

lmm.nointx <- lmer(mass.g ~ days.forelimb + water.level.reduc + (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)
lmm.nointx.log <- lmer(log(mass.g) ~ days.forelimb + water.level.reduc + (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)

lmm.null<- lmer(mass.g ~ (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)
lmm.null.log <- lmer(log(mass.g) ~ (1|clutch), data = filter(morph.data.mm, treatment == "low density"), na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, 
             lmm.nointx, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, 
             lmm.nointx.log,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T)
testDispersion(lmm.nointx)
testZeroInflation(lmm.nointx)
testCategorical(lmm.nointx, catPred = morph.data.mm$water.level.reduc[morph.data.mm$treatment == "low density"]) 

# check assumptions of best-fit model with log-transformed response variable
simulateResiduals(fittedModel = lmm.nointx.log, quantreg=T, plot = T)
testDispersion(lmm.nointx.log)
testZeroInflation(lmm.nointx.log)
testCategorical(lmm.nointx.log, catPred = morph.data.mm$water.level.reduc[morph.data.mm$treatment == "low density"]) 

final.mod = lmm.nointx

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)

#explanation of singular fit and why we can only assess for low density
View(morph.data.mm %>% 
       filter(treatment == "low density") %>%
       group_by(gs, gs.code, clutch, treatment, larv.tank.id, clutchtank, water.level.reduc) %>%
       summarise(n = n()))


# Supp Table 8a. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, and tank-level mean metamorphic mass on MASS ---------------------

# define candidate models
lmm.full <- lmer(mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log <- lmer(log(mass.g) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly2 <- lmer(mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly2 <- lmer(log(mass.g) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly3 <- lmer(mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly3 <- lmer(log(mass.g) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(mass.g ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(mass.g) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(mass.g ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(mass.g) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(mass.g ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(mass.g ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(mass.g) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log.poly3 <- lmer(log(mass.g) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(mass.g ~ 1 + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE) #singular fit
lmm.null.log<- lmer(log(mass.g) ~ 1 + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3, 
             lmm.nointxn.1, lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log,lmm.full.log.poly2, lmm.full.log.poly3,
                         lmm.nointxn.1.log, lmm.nointxn.log,lmm.nointxn.log.poly2,lmm.nointxn.log.poly3), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.poly3, form = morph.data.juv$treatment[is.na(morph.data.juv$mass.g)==FALSE])
plotResiduals(lmm.nointxn.poly3, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$mass.g)==FALSE])

simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment[is.na(morph.data.juv$mass.g)==FALSE])
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$mass.g)==FALSE])

final.mod = lmm.nointxn.log.poly3

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 8b. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, tank-level mean metamorphic mass, and tank-level mean body condition on BODY CONDITION (SCALED MASS INDEX, SMI) ---------------------

# define candidate models
lmm.full <- lmer(smi ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g*mean.mm.smi + scale(post.mm.weeks.num) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE) #singular fit with more complex and simplified random effect structure
lmm.full.log <- lmer(log(smi) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure

lmm.full.poly2 <- lmer(smi ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure
lmm.full.log.poly2 <- lmer(log(smi) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure

lmm.full.poly3 <- lmer(smi ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure
lmm.full.log.poly3 <- lmer(log(smi) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure

lmm.nointxn.1 <- lmer(smi ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g*mean.mm.smi + scale(post.mm.weeks.num) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(smi) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure

lmm.nointxn.1.poly2 <- lmer(smi ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g*mean.mm.smi + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.poly3 <- lmer(smi ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g*mean.mm.smi + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1.log.poly2 <- lmer(log(smi) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g*mean.mm.smi + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log.poly3 <- lmer(log(smi) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g*mean.mm.smi + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(smi ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + mean.mm.smi + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(smi) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure

lmm.nointxn.poly2 <- lmer(smi ~ treatment + mean.days.forelimb + mean.mm.mass.g + mean.mm.smi + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(smi ~ treatment + mean.days.forelimb + mean.mm.mass.g + mean.mm.smi + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(smi) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutch), data = morph.data.juv, na.action = na.omit, REML=FALSE)  #singular fit with more complex and simplified random effect structure
lmm.nointxn.log.poly3 <- lmer(log(smi) ~ treatment + mean.days.forelimb + mean.mm.mass.g + mean.mm.smi + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(smi ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.null.log<- lmer(log(smi) ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3,
             lmm.nointxn.1, lmm.nointxn.1.poly2, lmm.nointxn.1.poly3,
             lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, 
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log, lmm.full.log.poly2, lmm.full.log.poly3,
             lmm.nointxn.1.log, lmm.nointxn.1.log.poly2, lmm.nointxn.1.log.poly3,
             lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn.poly3 , quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.poly3 , form = morph.data.juv$treatment[is.na(morph.data.juv$mass.g)==FALSE])
plotResiduals(lmm.nointxn.poly3 , form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$mass.g)==FALSE])

simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment[is.na(morph.data.juv$mass.g)==FALSE])
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$mass.g)==FALSE])

final.mod = lmm.nointxn.poly3

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 8c. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, and tank-level mean metamorphic mass on SNOUT-VENT LENGTH ---------------------

# define candidate models
lmm.full <- lmer(svl.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log <- lmer(log(svl.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly2 <- lmer(svl.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly2 <- lmer(log(svl.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly3 <- lmer(svl.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly3 <- lmer(log(svl.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(svl.mm ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(svl.mm) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(svl.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(svl.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(svl.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(svl.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(svl.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log.poly3 <- lmer(log(svl.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE) 
lmm.null.log<- lmer(log(svl.mm) ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)  

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3, 
             lmm.nointxn.1, lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3,
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log,lmm.full.log.poly2, lmm.full.log.poly3,
             lmm.nointxn.1.log, lmm.nointxn.log,lmm.nointxn.log.poly2,lmm.nointxn.log.poly3,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

simulateResiduals(fittedModel = lmm.nointxn.log.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

final.mod = lmm.nointxn.log.poly2

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 8d. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, and tank-level mean metamorphic mass on FOREARM LENGTH ---------------------

# define candidate models
lmm.full <- lmer(r.forelimb.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log <- lmer(log(r.forelimb.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly2 <- lmer(r.forelimb.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly2 <- lmer(log(r.forelimb.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly3 <- lmer(r.forelimb.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly3 <- lmer(log(r.forelimb.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.forelimb.mm ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(r.forelimb.mm) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.forelimb.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(r.forelimb.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.forelimb.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(r.forelimb.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.forelimb.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log.poly3 <- lmer(log(r.forelimb.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.forelimb.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE) #singular fit
lmm.null.log<- lmer(log(r.forelimb.mm) ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3, 
             lmm.nointxn.1, lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log,lmm.full.log.poly2, lmm.full.log.poly3,
             lmm.nointxn.1.log, lmm.nointxn.log,lmm.nointxn.log.poly2,lmm.nointxn.log.poly3,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

simulateResiduals(fittedModel = lmm.nointxn.log.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

final.mod = lmm.nointxn.log.poly2

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 8e. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, and tank-level mean metamorphic mass on RIGHT TIBIA LENGTH ---------------------

# define candidate models
lmm.full <- lmer(r.tibia.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log <- lmer(log(r.tibia.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly2 <- lmer(r.tibia.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly2 <- lmer(log(r.tibia.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly3 <- lmer(r.tibia.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly3 <- lmer(log(r.tibia.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.tibia.mm ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(r.tibia.mm) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.tibia.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(r.tibia.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.tibia.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(r.tibia.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.tibia.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log.poly3 <- lmer(log(r.tibia.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.tibia.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.null.log<- lmer(log(r.tibia.mm) ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3, 
             lmm.nointxn.1, lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3,
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log,lmm.full.log.poly2, lmm.full.log.poly3,
             lmm.nointxn.1.log, lmm.nointxn.log,lmm.nointxn.log.poly2,lmm.nointxn.log.poly3,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

simulateResiduals(fittedModel = lmm.nointxn.log.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

final.mod = lmm.nointxn.log.poly2

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 8f. ANALYZE DATA - JUVENILE MORPHOMETRICS: Effect of rearing density, larval duration, and tank-level mean metamorphic mass on RIGHT THIGH LENGTH ---------------------

# define candidate models
lmm.full <- lmer(r.thigh.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log <- lmer(log(r.thigh.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly2 <- lmer(r.thigh.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly2 <- lmer(log(r.thigh.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.full.poly3 <- lmer(r.thigh.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.full.log.poly3 <- lmer(log(r.thigh.mm) ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.thigh.mm ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.1.log <- lmer(log(r.thigh.mm) ~ treatment + mean.days.forelimb*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.thigh.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log <- lmer(log(r.thigh.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.thigh.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.poly3 <- lmer(r.thigh.mm ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.thigh.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.nointxn.log.poly3 <- lmer(log(r.thigh.mm) ~ treatment + mean.days.forelimb + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.thigh.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)
lmm.null.log<- lmer(log(r.thigh.mm) ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# non-nested model selection using AICc for non-transformed response variable
arrange(AICc(lmm.full, lmm.full.poly2, lmm.full.poly3, 
             lmm.nointxn.1, lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3,
             lmm.null), AICc)

# non-nested model selection using AICc for log-transformed response variable
arrange(AICc(lmm.full.log,lmm.full.log.poly2, lmm.full.log.poly3,
             lmm.nointxn.1.log, lmm.nointxn.log,lmm.nointxn.log.poly2,lmm.nointxn.log.poly3,
             lmm.null.log), AICc)

# check assumptions of best-fit model with non-transformed response variable (lmm.nointxn.poly3) & log-transformed response variable (lmm.nointxn.log.poly3)
simulateResiduals(fittedModel = lmm.nointxn.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.poly3, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])


simulateResiduals(fittedModel = lmm.nointxn.log.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$treatment[is.na(morph.data.juv$treatment)==FALSE])
plotResiduals(lmm.nointxn.log.poly2, form = morph.data.juv$post.mm.weeks.num[is.na(morph.data.juv$treatment)==FALSE])

final.mod = lmm.nointxn.log.poly2

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# Supp Table 9a. ANALYZE DATA - TADPOLE SURVIVORSHIP: Effect of rearing density and age on % seeded tadpoles that survived weekly from start of experiment until metamorphosis ---------------------

# define candidate models
glmm.full.01 <- glmer(status ~ treatment*week + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.full.01.slopes <- glmer(status ~ treatment*week + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + week + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) #singular fit with more complex random effect, so simplified

glmm.null <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc 
model.sel = arrange(AICc(glmm.full.01, glmm.full.01.slopes,
             glmm.nointxn.01, glmm.nointxn.01.slopes,
             glmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(final.mod, form = survi.data.tad.exp.longform$treatment)
plotResiduals(final.mod, form = survi.data.tad.exp.longform$week)
plotResiduals(final.mod, form = survi.data.tad.exp.longform$water.level.reduc)
testDispersion(final.mod)

# final model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type III because interaction term significant
summary(final.mod)


# Supp Table 9b. ANALYZE DATA - JUVENILE SURVIVORSHIP: Effect of rearing density, age, and tank-level mean metamorphic mass on % seeded juveniles that survived weekly from metamorphosis until 8 months after metamorphosis ---------------------

# define candidate models
glmm.full.01 <- glmer(status ~ treatment*scale(mean.days.forelimb)*scale(mean.mm.mass.g) + postmm.week + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.full.01.slopes <- glmer(status ~ treatment*scale(mean.days.forelimb)*scale(mean.mm.mass.g) + postmm.week + (postmm.week||clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.1 <- glmer(status ~ treatment*scale(mean.days.forelimb) + mean.mm.mass.g + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.nointxn.01.1.slopes <- glmer(status ~ treatment*scale(mean.days.forelimb) + mean.mm.mass.g + postmm.week +  (postmm.week||clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.2 <- glmer(status ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.nointxn.01.2.slopes <- glmer(status ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + postmm.week +  (postmm.week||clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.3 <- glmer(status ~ treatment*mean.mm.mass.g + scale(mean.days.forelimb) + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.nointxn.01.3.slopes <- glmer(status ~ treatment*mean.mm.mass.g + scale(mean.days.forelimb) + postmm.week +  (postmm.week||clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + postmm.week + mean.days.forelimb + mean.mm.mass.g + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
glmm.nointxn.01.slopes <- glmer(status ~ treatment + postmm.week + mean.days.forelimb + mean.mm.mass.g + (postmm.week||clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc 
model.sel = arrange(AICc(glmm.full.01,glmm.full.01.slopes,
             glmm.nointxn.01.1, glmm.nointxn.01.1.slopes,
             glmm.nointxn.01.2, glmm.nointxn.01.2.slopes,
             glmm.nointxn.01.3, glmm.nointxn.01.3.slopes,
             glmm.nointxn.01, glmm.nointxn.01.slopes,
             glmm.null), AICc)

final.mod = eval(parse(text = paste(rownames(model.sel)[1])))

# check assumptions
simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(final.mod, form = survi.data.juv.longform$treatment)
plotResiduals(final.mod, form = survi.data.juv.longform$week)
testDispersion(final.mod)

# final model
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant
summary(final.mod)
