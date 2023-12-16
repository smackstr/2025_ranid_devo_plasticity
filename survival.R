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
library(DHARMa)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for survival data and development data
survi.data.tad = read.csv("Database_Survivorship - Tadpole Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.exp = read.csv("Database_Survivorship - Tadpole Experiment Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.juv = read.csv("Database_Survivorship - Froglet_Toadlet Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database
devo.data.nonmm = read.csv("Database_Metamorphosis - Non-Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database
morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
survi.data.tad = survi.data.tad[survi.data.tad$gs.code == "RS",]
survi.data.exp = survi.data.exp[survi.data.exp$gs.code == "RS",]
survi.data.juv = survi.data.juv[survi.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS" & is.na(devo.data$gs.code) == FALSE,]
#devo.data.nonmm = devo.data.nonmm[devo.data.nonmm$gs.code == "RS" & is.na(devo.data.nonmm$min.gosner.stage) == FALSE,]
morph.data.mm = morph.data.mm[morph.data.mm$gs.code == "RS",]

# subset databases to only include weeks up to 18
survi.data.juv = survi.data.juv[survi.data.juv$postmm.week <= 18 & is.na(survi.data.juv$postmm.week) == FALSE,]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
survi.data.tad$treatment[survi.data.tad$treatment == "control"] = "low density"
survi.data.exp$treatment[survi.data.exp$treatment == "control"] = "low density"
survi.data.juv$treatment[survi.data.juv$treatment == "control"] = "low density"
devo.data$treatment[devo.data$treatment == "control"] = "low density"
devo.data.nonmm$treatment[devo.data.nonmm$treatment == "control"] = "low density"
morph.data.mm$treatment[morph.data.mm$treatment == "control"] = "low density"

# subset devo dataset to only include first six individuals and non-overflow individuals
devo.data = devo.data[devo.data$first.six == "yes",]
devo.data = devo.data[devo.data$treatment != "overflow",]
morph.data.mm = morph.data.mm[morph.data.mm$first.six == "yes",]
morph.data.mm = morph.data.mm[morph.data.mm$treatment != "overflow",]

#change column classes
survi.data.tad$treatment = factor(survi.data.tad$treatment, levels = c("low density", "high density"))
survi.data.tad$larv.tank.id = factor(survi.data.tad$larv.tank.id)
survi.data.tad$water.level.reduc = factor(survi.data.tad$water.level.reduc)

survi.data.exp$treatment = factor(survi.data.exp$treatment, levels = c("low density", "high density"))
survi.data.exp$larv.tank.id = factor(survi.data.exp$larv.tank.id)
survi.data.exp$water.level.reduc = factor(survi.data.exp$water.level.reduc)

survi.data.juv$treatment = factor(survi.data.juv$treatment, levels = c("low density", "high density"))
survi.data.juv$juv.tank.id = factor(survi.data.juv$juv.tank.id)

devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment, levels = c("low density", "high density"))
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)
devo.data$days.forelimb = as.integer(devo.data$days.forelimb)
devo.data$days.forelimb.tail = as.integer(devo.data$days.forelimb.tail)

devo.data.nonmm$treatment = factor(devo.data.nonmm$treatment, levels = c("low density", "high density"))
devo.data.nonmm$larv.tank.id = factor(devo.data.nonmm$larv.tank.id)

morph.data.mm$first.six = factor(morph.data.mm$first.six, levels = c("yes", "no"))
morph.data.mm$treatment = factor(morph.data.mm$treatment, levels = c("low density", "high density"))
morph.data.mm$larv.tank.id = factor(morph.data.mm$larv.tank.id)
morph.data.mm$juv.tank.id = factor(morph.data.mm$juv.tank.id)
morph.data.mm$post.mm.weeks = factor(morph.data.mm$post.mm.weeks, ordered = TRUE, levels = c("0")) #set weeks as ordered factor

# create unique tank id for survi.data.tad, survi.data.exp, survi.data.juv, and devo.data
survi.data.tad$unique.id = paste(survi.data.tad$gs.code, survi.data.tad$clutch, survi.data.tad$larv.tank.id, sep = "_")
survi.data.exp$unique.id = paste(survi.data.exp$gs.code, survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_")
survi.data.juv$unique.id = paste(survi.data.juv$gs.code, survi.data.juv$clutch, survi.data.juv$juv.tank.id, sep = "_")
survi.data.juv$unique.id.wk = paste(survi.data.juv$gs.code, survi.data.juv$clutch, survi.data.juv$juv.tank.id, survi.data.juv$postmm.week, sep = "_")
devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.juv = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")
devo.data.nonmm$unique.id = paste(devo.data.nonmm$gs.code, devo.data.nonmm$clutch, devo.data.nonmm$larv.tank.id, sep = "_")
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")

#remove any NAs that have occurred
devo.data = devo.data[is.na(devo.data$gs) == FALSE,]

#create column to calculate average gosner stage
devo.data.nonmm$avg.gosner.stage = (devo.data.nonmm$min.gosner.stage + devo.data.nonmm$max.gosner.stage)/2

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


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to survi.data.exp -----------------------

# create column to store the developmental data, mean.days.forelimb represents days.forelimb for each individual (rather than an average, since we actually have individual level data for this) but need to keep column name the same as in the juvenile dataset so we can combine them in the next step
survi.data.exp$mean.days.forelimb = NA

# fill developmental data column for survi.data.tad.exp
for(i in 1:length(unique(survi.data.exp$unique.id))){
  survi.data.exp$mean.days.forelimb[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(survi.data.exp$unique.id)[i]], na.rm = TRUE)
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
  
  if(survi.data.juv$mean.days.forelimb[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] <= mean(devo.data$days.forelimb, na.rm = TRUE) - (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    survi.data.juv$devo.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "early"
  }
  
  if(survi.data.juv$mean.days.forelimb[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] >= mean(devo.data$days.forelimb, na.rm = TRUE) + (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    survi.data.juv$devo.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "late"
  }
}
survi.data.juv$devo.cat[is.na(survi.data.juv$devo.cat) == TRUE] = "mid"

#remove 0 weeks from survi.data.juv since we will only be analyzing survivorship after full weeks
survi.data.juv = survi.data.juv[survi.data.juv$postmm.week > 0,]


# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to survi.data.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
survi.data.juv$mean.mm.mass.g = NA

# fill developmental data column for survi.data.juv
for(i in 1:length(unique(survi.data.juv$unique.id))){
  survi.data.juv$mean.mm.mass.g[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(survi.data.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical size data (i.e. small, medium, large tanks)
survi.data.juv$mass.cat = NA

for(i in 1:length(unique(survi.data.juv$unique.id))){
  
  if(survi.data.juv$mean.mm.mass.g[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] <= mean(morph.data.mm$mass.g, na.rm = TRUE) - (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    survi.data.juv$mass.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "small"
  }
  
  if(survi.data.juv$mean.mm.mass.g[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]][1] >= mean(morph.data.mm$mass.g, na.rm = TRUE) + (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    survi.data.juv$mass.cat[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] = "large"
  }
}
survi.data.juv$mass.cat[is.na(survi.data.juv$mass.cat) == TRUE] = "med"


# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to survi.data.exp -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
survi.data.exp$mean.mm.mass.g = NA

#create unique id based on larv tank for morph.data.mm
morph.data.mm$unique.id.larv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, sep = "_")

# fill mean.mass.g for survi.data.exp -- currently calculating mean based on first six mass
for(i in 1:length(unique(survi.data.exp$unique.id))){
  survi.data.exp$mean.mm.mass.g[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id.larv ==  unique(survi.data.exp$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Scaled Mass Index (mean smi at mm) to survi.data.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
survi.data.juv$mean.mm.smi = NA

# fill developmental data column for survi.data.juv
for(i in 1:length(unique(survi.data.juv$unique.id))){
  survi.data.juv$mean.mm.smi[survi.data.juv$unique.id == unique(survi.data.juv$unique.id)[i]] <- mean(morph.data.mm$smi[morph.data.mm$unique.id ==  unique(survi.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATA: Generate survivorship dataframe in terms of 0 (dead) and 1 (alive) for tadpole weeks up to end of experiment ---------------

#create unique id that also has week attached to it
survi.data.tad.exp$unique.id.week = paste(survi.data.tad.exp$unique.id, survi.data.tad.exp$week, sep = "_")

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.tad.exp[survi.data.tad.exp$unique.id.week == unique(survi.data.tad.exp$unique.id.week)[1],]
temp3 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                   gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                   gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                   clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                   water.level.reduc = rep(temp$water.level.reduc, temp$seed.num - temp$leth.samp.num.cumul),
                   status = NA
)
temp3$status[0:(temp$photo.num + temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id.week))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id.week == unique(survi.data.tad.exp$unique.id.week)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                     gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                     water.level.reduc = rep(temp$water.level.reduc, temp$seed.num - temp$leth.samp.num.cumul),
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
                   mean.mm.mass.g = rep(temp$mean.mm.mass.g, temp$seed.num - temp$leth.samp.num),
                   mean.mm.smi = rep(temp$mean.mm.smi, temp$seed.num - temp$leth.samp.num),
                   devo.cat = rep(temp$devo.cat, temp$seed.num - temp$leth.samp.num),
                   unique.id = rep(temp$unique.id, temp$seed.num - temp$leth.samp.num),
                   status = NA
)
temp3$status[0:(temp$live.num)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
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
                     mean.mm.mass.g = rep(temp$mean.mm.mass.g, temp$seed.num - temp$leth.samp.num),
                     mean.mm.smi = rep(temp$mean.mm.smi, temp$seed.num - temp$leth.samp.num),
                     devo.cat = rep(temp$devo.cat, temp$seed.num - temp$leth.samp.num),
                     unique.id = rep(temp$unique.id, temp$seed.num - temp$leth.samp.num),
                     status = NA
  )
  temp2$status[0:(temp$live.num)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
survi.data.juv.longform = temp3
rm(temp3)


# COMPILE DATASETS: Create clutchtank column we can use for random effect -----------------------
devo.data$clutchtank = factor(paste(devo.data$clutch, devo.data$larv.tank.id, sep = "_"))
survi.data.exp$clutchtank = factor(paste(survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_"))
survi.data.tad$clutchtank = factor(paste(survi.data.tad$clutch, survi.data.tad$larv.tank.id, sep = "_"))
survi.data.tad.exp$clutchtank = factor(paste(survi.data.tad.exp$clutch, survi.data.tad.exp$larv.tank.id, sep = "_"))
survi.data.juv$clutchtank = factor(paste(survi.data.juv$clutch, survi.data.juv$juv.tank.id, sep = "_"))
survi.data.tad.exp.longform$clutchtank = factor(paste(survi.data.tad.exp.longform$clutch, survi.data.tad.exp.longform$larv.tank.id, sep = "_"))
survi.data.juv.longform$clutchtank = factor(paste(survi.data.juv.longform$clutch, survi.data.juv.longform$larv.tank.id, sep = "_"))


# ANALYZE DATA: Effect of rearing density, water level reduction, week, and larv tank id nested within clutch on % surviving individuals weekly up until metamorphosis ---------------------
# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but maximum likelihood does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.

# model definition 
lmm.full <- lmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, weights = (seed.num-leth.samp.num.cumul))

glmm.full <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.probit <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial(link="probit"), weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.slopes <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn.slopes <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.null <- glmer((photo.num + metamorph.num.cumul)/(seed.num-leth.samp.num.cumul) ~ (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))


# check assumptions
simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.probit, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.slopes, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.nointxn, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.nointxn.slopes, quantreg=T, plot = T) #provides summary of model fitting tests

#tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
plotResiduals(lmm.full, form = survi.data.tad.exp$treatment)
plotResiduals(glmm.full, form = survi.data.tad.exp$treatment)
plotResiduals(glmm.full.probit, form = survi.data.tad.exp$treatment)
plotResiduals(glmm.nointxn, form = survi.data.tad.exp$treatment)
plotResiduals(glmm.nointxn.slopes, form = survi.data.tad.exp$treatment)

plotResiduals(lmm.full, form = survi.data.tad.exp$week)
plotResiduals(glmm.full, form = survi.data.tad.exp$week)
plotResiduals(glmm.full.probit, form = survi.data.tad.exp$week)
plotResiduals(glmm.nointxn, form = survi.data.tad.exp$week)
plotResiduals(glmm.nointxn.slopes, form = survi.data.tad.exp$week)

plotResiduals(lmm.full, form = survi.data.tad.exp$water.level.reduc)
plotResiduals(glmm.full, form = survi.data.tad.exp$water.level.reduc)
plotResiduals(glmm.full.slopes, form = survi.data.tad.exp$water.level.reduc)
plotResiduals(glmm.full.probit, form = survi.data.tad.exp$water.level.reduc)
plotResiduals(glmm.nointxn, form = survi.data.tad.exp$water.level.reduc)
plotResiduals(glmm.nointxn.slopes, form = survi.data.tad.exp$water.level.reduc)

testDispersion(lmm.full)
testDispersion(glmm.full) #tests for over- and under-dispersion
testDispersion(glmm.full.slopes) #tests for over- and under-dispersion
testDispersion(glmm.full.probit) #tests for over- and under-dispersion
testDispersion(glmm.nointxn) #tests for over- and under-dispersion
testDispersion(glmm.nointxn.slopes) #tests for over- and under-dispersion

#examine the plots for for each tank, treatment, and clutch, regress survival proportion by week using linear and then using linear and you can see how much beter the glm is and how it does look like distinct slopes for high density vs low density
ggplot(survi.data.tad.exp, 
       aes(x = week, y = prop.survi, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = survi.data.tad.exp, method = "lm", aes(x=week, y=prop.survi), inherit.aes = F, se = F, color="black") +
  geom_smooth(data = survi.data.tad.exp, method = "glm", method.args = list(family = "binomial"), aes(x=week, y=prop.survi), inherit.aes = F, se = F, color="black")

ggplot() +
  stat_summary(data = survi.data.tad.exp, fun.y=mean, geom="line", size = 0.5, alpha = 0.5, aes(x = week, y = prop.survi, color = treatment, group = clutchtank)) +
  stat_summary(data = survi.data.tad.exp, fun.y=mean, geom="point", color = "black", pch=21, alpha = 0.8, size=2, aes(x = week, y = prop.survi, fill=treatment)) +
  #geom_smooth(data = survi.data.tad.exp, method = "lm", aes(x=week, y=prop.seed.forelimb), inherit.aes = F, se = F, color="black")
  geom_smooth(data = survi.data.tad.exp, method = "glm", method.args = list(family = "binomial"(link = "logit")), aes(x=week, y=prop.survi, color = treatment, group = treatment, linetype = treatment), inherit.aes = F, se = F, size = 1.5) +
  scale_x_continuous(limits = c(1,11), breaks = seq(1,11,1))

#model selection using AIC among the logit link models. assess of modeling week as polynomial is better fit within binomial glmm - yes some support for this
anova(glmm.full, glmm.full.poly2, glmm.full.poly3, glmm.full.slopes, glmm.full.slopes.poly2, glmm.full.slopes.poly3, glmm.nointxn, glmm.nointxn.slopes)

#but assumption plots indicate that we have a problem with overdispersion, which occurs when error (residuals) are more variable than expected from the theorized distribution, so modeling (0,1) data to see if that helps. Note some models estimate intercept and slope, separately, by random factor: (1 | random.factor) + (0 + fixed.factor | random.factor). An alternative way to write this is using the double-bar notation fixed.factor + (fixed.factor || random.factor). when graphing the data, while it looks like a week^2 or week^3 could be helpful, this is likely just because we lose clutches b and c for week 11 so the mean is artificially brought up of this, keeping the formula linear.
glmm.full.01 <- glmer(status ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes <- glmer(status ~ treatment*week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null <- glmer(status ~ (1|clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.01.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.01.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.01.slopes, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.01.slopes.poly2, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full.01.slopes.poly3, quantreg=T, plot = T) #provides summary of model fitting tests

#model selection out of 0/1 logit models using likelihood-ratio test
anova(glmm.full.01, glmm.nointxn.01, glmm.nointxn.01.slopes, glmm.full.01.slopes, glmm.null)
anova(glmm.full.01, glmm.nointxn.01, glmm.null)

#model selection out of best-fit 0/1 logit models using AICc
AICc(glmm.nointxn.01, glmm.nointxn.01.slopes)

#doublecheck assumptions
simulateResiduals(fittedModel = glmm.nointxn.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.nointxn.01, form = survi.data.tad.exp.longform$treatment)
plotResiduals(glmm.nointxn.01, form = survi.data.tad.exp.longform$week)
plotResiduals(glmm.nointxn.01, form = survi.data.tad.exp.longform$water.level.reduc)
testDispersion(glmm.nointxn.01)

# Final Model:
final.mod = glmm.nointxn.01
Anova(final.mod)
summary(final.mod)

# create dataframe of predicted values that can be plotted on ggplot later
predicted.df.prop.survi <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, larval duration, week, mean mass at metamorphosis, and juv tank id nested within clutch on % surviving individuals weekly AFTER metamorphhosis---------------------

#examine the plots for for each tank, treatment, and clutch, regress survival proportion by week using linear and then using linear and you can see how much beter the glm is and how it does look like distinct slopes for high density vs low density
ggplot(survi.data.juv, 
       aes(x = postmm.week, y = prop.survi, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = survi.data.juv, method = "lm", aes(x=postmm.week, y=prop.survi), inherit.aes = F, se = F, color="black") +
  geom_smooth(data = survi.data.juv, method = "glm", method.args = list(family = "binomial"), aes(x=postmm.week, y=prop.survi), inherit.aes = F, se = F, color="black")

ggplot() +
  stat_summary(data = survi.data.juv, fun.y=mean, geom="line", size = 0.5, alpha = 0.5, aes(x = postmm.week, y = prop.survi, color = treatment, group = clutchtank)) +
  stat_summary(data = survi.data.juv, fun.y=mean, geom="point", color = "black", pch=21, alpha = 0.8, size=2, aes(x = postmm.week, y = prop.survi, fill=treatment)) +
  #geom_smooth(data = survi.data.tad.exp, method = "lm", aes(x=week, y=prop.seed.forelimb), inherit.aes = F, se = F, color="black")
  geom_smooth(data = survi.data.juv, method = "glm", method.args = list(family = "binomial"(link = "logit")), aes(x=postmm.week, y=prop.survi, color = treatment, group = treatment, linetype = treatment), inherit.aes = F, se = F, size = 1.5) +
  scale_x_continuous(limits = c(1,11), breaks = seq(1,11,1))

# model definition - setting one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. since experimental study, not testing the inclusion or exclusion of fixed effects, but rather assessing what interactions among fixed effects to include in final model.Note some models estimate intercept and slope, separately, by random factor: (1 | random.factor) + (0 + fixed.factor | random.factor). An alternative way to write this is using the double-bar notation fixed.factor + (fixed.factor || random.factor).
glmm.full.01 <- glmer(status ~ treatment*scale(mean.days.forelimb)*scale(mean.mm.mass.g) + postmm.week + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.2 <- glmer(status ~ treatment*scale(mean.days.forelimb)*scale(mean.mm.mass.g)*mean.mm.smi + postmm.week + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.1 <- glmer(status ~ treatment*scale(mean.days.forelimb) + mean.mm.mass.g + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.2 <- glmer(status ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.3 <- glmer(status ~ treatment*mean.mm.mass.g + scale(mean.days.forelimb) + postmm.week +  (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + postmm.week + scale(mean.days.forelimb) + mean.mm.mass.g + (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null <- glmer(status ~ (1|clutchtank), data=survi.data.juv.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

#model selection out of 0/1 logit models
anova(glmm.full.01, glmm.full.01.2, glmm.nointxn.01.1, glmm.nointxn.01.2, glmm.nointxn.01.3, glmm.nointxn.01, glmm.null)
anova(glmm.full.01, glmm.full.01.sigmoid, glmm.null)

#doublecheck assumptions
simulationOutput <- simulateResiduals(fittedModel = glmm.full.01.2, quantreg=T, plot = T) #provides summary of model fitting tests
testOutliers(simulationOutput, type = 'bootstrap') #recommended re-test with type=outliers indicates no issues
plotResiduals(simulationOutput, form = survi.data.juv.longform$treatment)
plotResiduals(simulationOutput, form = survi.data.juv.longform$postmm.week)
testDispersion(simulationOutput)

#SMS YOU ARE HERE: currently figuring out if it makes sense to even run glmm.full.01.2 since mean.days.forelimb, mean.mm.mass.g, and mean.mm.smi are all correlated


# Final Model:
final.mod = glmm.full.01.2
Anova(final.mod, type = "II")
summary(final.mod)
cov2cor(vcov(final.mod)) #assess correlation matrix between fixed effects

# create dataframe of predicted values that can be plotted on ggplot later
predicted.df.juv.prop.survi <- ggeffects::ggpredict(final.mod, terms = c("postmm.week[all]", "treatment"), type = "random", interval = "confidence")

# significant three-way interaction suggests mean days forelimb and mean mm size affect survivorship differently in high and low density so splitting dataset into two different analysis based on treatment to understand better
glmm.full.01.ld <- glmer(status ~ scale(mean.days.forelimb)*mean.mm.mass.g + postmm.week + (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.ld <- glmer(status ~ scale(mean.days.forelimb) + mean.mm.mass.g + postmm.week + (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.ld <- glmer(status ~ (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.hd <- glmer(status ~ scale(mean.days.forelimb)*mean.mm.mass.g + postmm.week + (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "high density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.hd <- glmer(status ~ scale(mean.days.forelimb) + mean.mm.mass.g + postmm.week + (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "high density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.hd <- glmer(status ~ (1|clutchtank), data=survi.data.juv.longform[survi.data.juv.longform$treatment == "high density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

anova(glmm.full.01.ld, glmm.nointxn.01.ld, glmm.null.ld)
anova(glmm.full.01.hd, glmm.nointxn.01.hd, glmm.null.hd)

final.mod.ld = glmm.nointxn.01.ld
final.mod.hd = glmm.full.01.hd

Anova(final.mod.ld, type = 'II')
Anova(final.mod.hd, type = 'II')

summary(final.mod.ld)
summary(final.mod.hd) #at larger sizes, mean days forelimb has less of a negative effect


# ANALYZE DATA: Effect of rearing density, larval duration, week, and juv tank id nested within clutch on % surviving individuals weekly BUT NOW USING CATEGORY FOR SIZE ---------------------

# model definition - setting one random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. Subsetting weeks 1-12 to represent breadth of data we are collecting, but will likely change this once we get another set of morphometric results for weeks 14-16

#need this to be prop.survi so that predict function will work for graph
glmm.nointxn <- glmer(prop.survi ~ treatment + mass.cat + postmm.week + (1|clutch:juv.tank.id), data=survi.data.juv, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num))

# # create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.surv <- ggpredict(glmm.nointxn, terms = c("postmm.week", "mass.cat"), type = "random", interval = "confidence") #ggemmeans returns confidence intervals while ggpredict does not...suspect this has something to do with categorical variable combined with logistic regression? but tutorials seem to be doing the exact same thing I'm doing and get confidence intervals. Could also be because expected 0 or 1 but getting a proportion?


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
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
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
  scale_y_continuous(name = "percent survived", limits = c(85, 100)) +
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

# percent survived from tank seeding to tank close-out without clutch included but now including predicted fit from model and plotting all tanks across time
ggplot() + 
  geom_jitter(width = 0.4, size = 2, alpha = 0.7, data = survi.data.tad.exp, pch = 21, aes(y=prop.survi, x = week, color = treatment, fill = treatment)) +
  geom_ribbon(data = predicted.df.prop.survi,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, fill = group), show.legend = FALSE) +
  geom_line(data = predicted.df.prop.survi,
            alpha = 1,
            size = 1.8,
            mapping = aes(x = x, y = predicted, color = group), show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent survivorship", limits=c(0,1.0), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), labels = c(0,10,20,30,40,50,60,70,80,90,100)) +
  scale_x_continuous(name = "pre-metamorphic age (weeks)", breaks = seq(0, 12, by = 1))


# PLOT DATASETS: Effect of predictors on survival after metamorphosis ---------------------

# x-y plot with smoothed binomial curve WITHOUT clutch
masseffect <- ggplot() + 
  facet_grid(cols=vars(treatment)) + 
  geom_point(size = 1, alpha = 0.7, data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), group=clutchtank, color=mean.mm.mass.g)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "glm", method.args = list(family = "binomial"), data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), color = mean.mm.mass.g, group=clutchtank)) +
  # geom_line(data = predicted.df.juv.prop.survi,
  #           alpha = 1,
  #           size = 1.8,
  #           mapping = aes(x = x, y = predicted, linetype = group), color = "darkred", show.legend = TRUE) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_fill_manual(values = c("gray1", "gray40", "gray85")) +
  scale_linetype_manual(values=c(1,3)) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent juvenile survivorship", labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "post-metamorphic age (week)", breaks = seq(1, 18, by = 1))

dayseffect <- ggplot() + 
  facet_grid(cols=vars(treatment)) + 
  geom_point(size = 1, alpha = 0.7, data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), group=clutchtank, color=mean.days.forelimb)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "glm", method.args = list(family = "binomial"), data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), color = mean.days.forelimb, group=clutchtank)) +
  # geom_line(data = predicted.df.juv.prop.survi,
  #           alpha = 1,
  #           size = 1.8,
  #           mapping = aes(x = x, y = predicted, linetype = group), color = "darkred", show.legend = TRUE) +
  scale_color_gradient(low = "gray1", high = "gray85") +
  scale_fill_manual(values = c("gray85", "gray40", "gray1")) +
  scale_linetype_manual(values=c(1,3)) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent juvenile survivorship", labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "post-metamorphic age (week)", breaks = seq(1, 18, by = 1))


ggplot() + 
  facet_grid(cols=vars(treatment)) + 
  geom_point(size = 1, alpha = 0.7, data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), group=clutchtank, color=mean.mm.mass.g)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "glm", method.args = list(family = "binomial"), data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), color = mean.mm.mass.g, group=clutchtank)) +
  geom_line(data = predicted.df.juv.prop.survi,
            alpha = 1,
            size = 1.8,
            mapping = aes(x = x, y = predicted, linetype = group), color = "darkred", show.legend = TRUE) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_fill_manual(values = c("gray1", "gray40", "gray85")) +
  
  new_scale_color() +
  new_scale_fill() +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
scale_fill_manual(values = c("darkblue", "cornflowerblue", "lightblue")) +
  geom_point(size = 1, alpha = 0.7, data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), group=clutchtank, color=mean.days.forelimb)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "glm", method.args = list(family = "binomial"), data = survi.data.juv, aes(y=prop.survi, x = as.numeric(postmm.week), color = mean.days.forelimb, group=clutchtank)) +
  
  scale_linetype_manual(values=c(1,3)) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "percent juvenile survivorship", labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "post-metamorphic age (week)", breaks = seq(1, 18, by = 1))

ggarrange(masseffect, dayseffect,
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))



# PLOT DATA: Create panel plot with all survival data across all sampling points --------------
ggarrange(plot.survi1, plot.survi2, plot.survi3,
          ncol = 3,
          nrow = 1,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b", "c"),
          font.label = list(size = 20, color = "black"))


# PLOT DATA: Examine mean metamorphic mass and mean larval duration across rearing density to understand three-way interaction with juvenile survivorship across week --------------

ggplot(data = survi.data.juv[survi.data.juv$postmm.week == 1,]) +
  geom_point(aes(x=treatment, y=mean.mm.mass.g, color = treatment, group = clutchtank)) +
  stat_summary(fun = "mean", geom="point", color = "black", pch=21, size=4,
               aes(x=treatment, y=mean.mm.mass.g, fill = treatment)) +
  geom_boxplot(alpha = 0.7,
               aes(x=treatment, y=mean.mm.mass.g, color = treatment))

ggplot(data = survi.data.juv[survi.data.juv$postmm.week == 1,]) +
  geom_point(aes(x=treatment, y=mean.days.forelimb, color = treatment, group = clutchtank)) +
  stat_summary(fun = "mean", geom="point", color = "black", pch=21, size=4,
               aes(x=treatment, y=mean.days.forelimb, fill = treatment)) +
  geom_boxplot(alpha = 0.7,
               aes(x=treatment, y=mean.days.forelimb, color = treatment))


# PLOT DATA: Examine relationship between metamorphic mass and larval duration across rearing densityto understand three-way interaction with juvenile survivorship across week --------------

ggplot(data = morph.data.mm) +
  facet_grid(cols=vars(treatment)) +
  geom_jitter(aes(x=days.forelimb, y=mass.g, color = clutch)) +
  geom_smooth(se=F, method = "lm", aes(x=days.forelimb, y=mass.g, group = clutchtank, color = clutch))


ggplot(data = survi.data.juv[survi.data.juv$postmm.week == 1,]) +
  geom_point(aes(x=treatment, y=mean.days.forelimb, color = treatment, group = clutchtank)) +
  stat_summary(fun = "mean", geom="point", color = "black", pch=21, size=4,
               aes(x=treatment, y=mean.days.forelimb, fill = treatment)) +
  geom_boxplot(alpha = 0.7,
               aes(x=treatment, y=mean.days.forelimb, color = treatment))


