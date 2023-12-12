# this script creates figures of the developmental results for the 2023 ranid developmental paper
# includes effects of rearing density on developmental timepoints (1. time to forelimb emergence, time to complete metamorphosis, laterality of forelimb emergence)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)
library(DHARMa)

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

# subset databases to only include first six individuals and non-overflow
devo.data = devo.data[devo.data$treatment != "overflow",]
devo.data = devo.data[devo.data$first.six == "yes",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data$treatment[devo.data$treatment == "control"] = "low density"
survi.data.exp$treatment[survi.data.exp$treatment == "control"] = "low density"
survi.data.tad$treatment[survi.data.tad$treatment == "control"] = "low density"

#change column classes
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment, levels = c("low density", "high density"))
devo.data$clutch = factor(devo.data$clutch)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)
devo.data$days.forelimb = as.integer(devo.data$days.forelimb)
devo.data$days.forelimb.tail = as.integer(devo.data$days.forelimb.tail)

survi.data.exp$treatment = factor(survi.data.exp$treatment, levels = c("low density", "high density"))
survi.data.exp$larv.tank.id = factor(survi.data.exp$larv.tank.id)
survi.data.exp$prop.surv.forelimb = as.numeric(survi.data.exp$prop.surv.forelimb)
survi.data.exp$water.level.reduc = factor(survi.data.exp$water.level.reduc)

survi.data.tad$treatment = factor(survi.data.tad$treatment, levels = c("low density", "high density"))
survi.data.tad$larv.tank.id = factor(survi.data.tad$larv.tank.id)
survi.data.tad$prop.surv.forelimb = as.numeric(survi.data.tad$prop.surv.forelimb)
survi.data.tad$water.level.reduc = factor(survi.data.tad$water.level.reduc)

survi.data.tad$week = factor(survi.data.tad$week, ordered = TRUE, levels = seq(1, 10, by = 1)) # week as an ordinal factor
survi.data.exp$week = factor(survi.data.exp$week, ordered = TRUE, levels = c(10,11)) # week as an ordinal factor

# create unique tank id for survi.data.exp and devo.data
survi.data.exp$unique.id = paste(survi.data.exp$gs.code, survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_")
devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, devo.data$animal.id, sep = "_")
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


# COMPILE DATA: Generate end of experiment metamorphosis dataframe in terms of 0 (tadpole) and 1 (metamorphosed) -----------------------
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
morph.data.exp.longform = temp3
rm(temp3)


# COMPILE DATASETS: Add Developmental Category (devo.cat) to devo.data -----------------------

# create column to store categorical developmental data (i.e. early, mid, late developers)

devo.data$devo.cat = NA

devo.data$devo.cat[devo.data$days.forelimb <= mean(devo.data$days.forelimb, na.rm = TRUE) - (sd(devo.data$days.forelimb, na.rm = TRUE)/2)] = "early"
devo.data$devo.cat[devo.data$days.forelimb >= mean(devo.data$days.forelimb, na.rm = TRUE) + (sd(devo.data$days.forelimb, na.rm = TRUE)/2)] = "late"
devo.data$devo.cat[is.na(devo.data$devo.cat) == TRUE] = "mid"


# COMPILE DATA: Generate end of experiment metamorphosis dataframe in terms of 0 (tadpole) and 1 (metamorph) BUT NOW ONLY CONSIDERING SURVIVING INDIVIDUALS -----------------------
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
morph.data.exp.longformsurvi = temp3
rm(temp3)

# COMPILE DATASETS: Combine tadpole and end of experiment dataframes -----------------------
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


# COMPILE DATA: Generate weekly metamorphosis dataframe in terms of 0 (tadpole) and 1 (metamorphosed) -----------------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (metamorphosed) is calculated from the metamorphosed individuals
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
temp3$status[0:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id.week))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id.week == unique(survi.data.tad.exp$unique.id.week)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$seed.num - temp$leth.samp.num.cumul),
                    gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num),
                     water.level.reduc = rep(temp$water.level.reduc, temp$seed.num - temp$leth.samp.num),
                     status = NA
  )
  temp2$status[0:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
morph.data.tad.exp.longform = temp3
rm(temp3)
morph.data.tad.exp.longform$status = factor(morph.data.tad.exp.longform$status)
morph.data.tad.exp.longform$water.level.reduc = factor(morph.data.tad.exp.longform$water.level.reduc)

# COMPILE DATA: Generate weekly metamorphosis dataframe in terms of 0 (tadpole) and 1 (metamorphosed) BUT ONLY CONSIDERING SURVIVING INDIVIDUALS OUT OF SEEDED -----------------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.tad.exp[survi.data.tad.exp$unique.id.week == unique(survi.data.tad.exp$unique.id.week)[1],]
temp3 = data.frame(week = rep(temp$week, temp$photo.num + temp$metamorph.num.cumul),
                   gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                   gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                   clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                   treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                   water.level.reduc = rep(temp$water.level.reduc, temp$photo.num + temp$metamorph.num.cumul),
                   status = NA
)
temp3$status[0:(temp$metamorph.num.cumul)] = 1
temp3$status[is.na(temp3$status) == TRUE] = 0
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's status for each tank and for each week)
for(i in 2:length(unique(survi.data.tad.exp$unique.id.week))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id.week == unique(survi.data.tad.exp$unique.id.week)[i],]
  temp2 = data.frame(week = rep(temp$week, temp$photo.num + temp$metamorph.num.cumul),
                     gs = rep(temp$gs, temp$photo.num + temp$metamorph.num.cumul),
                     gs.code = rep(temp$gs.code, temp$photo.num + temp$metamorph.num.cumul),
                     clutch = rep(temp$clutch, temp$photo.num + temp$metamorph.num.cumul),
                     treatment = rep(temp$treatment, temp$photo.num + temp$metamorph.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$photo.num + temp$metamorph.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$photo.num + temp$metamorph.num.cumul),
                     water.level.reduc = rep(temp$water.level.reduc, temp$photo.num + temp$metamorph.num.cumul),
                     status = NA
  )
  temp2$status[0:(temp$metamorph.num.cumul)] = 1
  temp2$status[is.na(temp2$status) == TRUE] = 0
  temp3 = rbind(temp3, temp2)
  rm(temp,temp2)
}
morph.data.tad.exp.longformsurvi = temp3
rm(temp3)
morph.data.tad.exp.longformsurvi$status = factor(morph.data.tad.exp.longformsurvi$status)
morph.data.tad.exp.longformsurvi$water.level.reduc = factor(morph.data.tad.exp.longformsurvi$water.level.reduc)

# COMPILE DATASETS: Add Gosner staging (mean devo stage) to survi.data.exp -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
survi.data.exp$tank.devo.stage = NA

# fill developmental data column for survi.data.exp.
for(i in 1:length(unique(survi.data.exp$unique.id))){
  
  #create vector that includes all the ones that already metamorphosed as stage 42 + non-metamorphosed stage
  temp = c(rep(42, survi.data.exp$metamorph.num.cumul[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]]), devo.data.nonmm$avg.gosner.stage[devo.data.nonmm$unique.id ==  unique(survi.data.exp$unique.id)[i]])
  
  #find mean of vector. purposefully not adding na.rm = TRUE so temporarily subset out the tanks we haven't staged yet
  survi.data.exp$tank.devo.stage[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] = mean(temp)
}


# COMPILE DATASETS: Create individual-level dataset for Gosner staging so we can visualize the distribution -----------------------

#set up database to fill
#number of rows (=number of tadpoles) is calculated from the seed number MINUS the lethal sampled number
#number of 1 (alive) is calculated from the photo number PLUS metamorphosed individuals
temp = survi.data.exp[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[1],]
temp = temp[temp$week == max(temp$week),]
temp2 = devo.data.nonmm[devo.data.nonmm$unique.id == unique(survi.data.exp$unique.id)[1],]
temp3 = data.frame(gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                   gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                   clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                   treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                   larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                   devo.stage = NA
)
temp3$devo.stage[0:temp$metamorph.num.cumul] = 42
temp3$devo.stage[(temp$metamorph.num.cumul + 1):((temp$metamorph.num.cumul) + length(temp2$avg.gosner.stage))] = temp2$avg.gosner.stage
rm(temp)

#fill in database with remainder of unique ids (representing each tadpole's devo stage for each tank )
for(i in 2:length(unique(survi.data.exp$unique.id))){
  temp = survi.data.tad.exp[survi.data.tad.exp$unique.id == unique(survi.data.exp$unique.id)[i],]
  temp = temp[temp$week == max(temp$week),]
  temp2 = devo.data.nonmm[devo.data.nonmm$unique.id == unique(survi.data.exp$unique.id)[i],]
  temp4 = data.frame(gs = rep(temp$gs, temp$seed.num - temp$leth.samp.num.cumul),
                     gs.code = rep(temp$gs.code, temp$seed.num - temp$leth.samp.num.cumul),
                     clutch = rep(temp$clutch, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment = rep(temp$treatment, temp$seed.num - temp$leth.samp.num.cumul),
                     treatment.code = rep(temp$treatment.code, temp$seed.num - temp$leth.samp.num.cumul),
                     larv.tank.id = rep(temp$larv.tank.id, temp$seed.num - temp$leth.samp.num.cumul),
                     devo.stage = NA
  )
  temp4$devo.stage[0:temp$metamorph.num.cumul] = 42
  temp4$devo.stage[(temp$metamorph.num.cumul + 1):((temp$metamorph.num.cumul) + length(temp2$avg.gosner.stage))] = temp2$avg.gosner.stage
  temp3 = rbind(temp3, temp4)
  rm(temp,temp2,temp4)
}

data.gstage = temp3
rm(temp3)

# COMPILE DATASETS: Create clutchtank column we can use for random effect -----------------------
devo.data$clutchtank = factor(paste(devo.data$clutch, devo.data$larv.tank.id, sep = "_"))
survi.data.exp$clutchtank = factor(paste(survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_"))
survi.data.tad$clutchtank = factor(paste(survi.data.tad$clutch, survi.data.tad$larv.tank.id, sep = "_"))
survi.data.tad.exp$clutchtank = factor(paste(survi.data.tad.exp$clutch, survi.data.tad.exp$larv.tank.id, sep = "_"))
survi.data.tad.exp.longform$clutchtank = factor(paste(survi.data.tad.exp.longform$clutch, survi.data.tad.exp.longform$larv.tank.id, sep = "_"))
survi.data.tad.exp.longformsurvi$clutchtank = factor(paste(survi.data.tad.exp.longformsurvi$clutch, survi.data.tad.exp.longformsurvi$larv.tank.id, sep = "_"))
morph.data.tad.exp.longform$clutchtank = factor(paste(morph.data.tad.exp.longform$clutch, morph.data.tad.exp.longform$larv.tank.id, sep = "_"))
morph.data.tad.exp.longformsurvi$clutchtank = factor(paste(morph.data.tad.exp.longformsurvi$clutch, morph.data.tad.exp.longformsurvi$larv.tank.id, sep = "_"))


# COMPILE DATASETS: Create 10% dataset to identify the first 10% to metamorphose for each rearing density -----------------------
devo.data.first10 = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
devo.data.first10 = devo.data.first10[devo.data.first10$gs.code == "RS",]

# subset databases to only include non-overflow
devo.data.first10 = devo.data.first10[devo.data.first10$treatment != "overflow",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data.first10$treatment[devo.data.first10$treatment == "control"] = "low density"

#change column classes
devo.data.first10$first.six = factor(devo.data.first10$first.six, levels = c("yes", "no"))
devo.data.first10$treatment = factor(devo.data.first10$treatment, levels = c("low density", "high density"))
devo.data.first10$clutch = factor(devo.data.first10$clutch)
devo.data.first10$larv.tank.id = factor(devo.data.first10$larv.tank.id)
devo.data.first10$juv.tank.id = factor(devo.data.first10$juv.tank.id)
devo.data.first10$days.forelimb = as.integer(devo.data.first10$days.forelimb)
devo.data.first10$days.forelimb.tail = as.integer(devo.data.first10$days.forelimb.tail)

# create unique.id for devo.data.first10
devo.data.first10$unique.id = paste(devo.data.first10$gs.code, devo.data.first10$clutch, devo.data.first10$larv.tank.id, sep = "_")
devo.data.first10$unique.id.indiv = paste(devo.data.first10$gs.code, devo.data.first10$clutch, devo.data.first10$larv.tank.id, devo.data.first10$animal.id, sep =

devo.data.first10$first.10perc = "NA"
devo.data.first10$first.10perc[devo.data.first10$treatment == "high density" & devo.data.first10$animal.id <=10] = "yes"
devo.data.first10$first.10perc[devo.data.first10$treatment == "high density" & devo.data.first10$animal.id >10] = "no" #but need to not filter out by yes to do this
devo.data.first10$first.10perc[devo.data.first10$treatment == "low density" & devo.data.first10$animal.id <=6] = "yes"
devo.data.first10$first.10perc[devo.data.first10$treatment == "low density" & devo.data.first10$animal.id >6] = "no"

#how many high density tanks did we remove that hadn't reached 10?
first10.summ = devo.data.first10 %>%
  group_by(treatment, clutch, larv.tank.id) %>%
  summarize(n = n())
nrow(temp[first10.summ$treatment == "high density" & first10.summ$n < 10,]) #this provides number of tanks that we are removing by opportunistically looking at first 10% of tanks rather than first six individuals

#create clutchtank for later analysis
devo.data.first10$clutchtank = paste(devo.data.first10$clutch, devo.data.first10$larv.tank.id, sep="_")


# ANALYZE DATA: Effect of rearing density and clutch:larval tank on larval duration for FIRST SIX ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but maximum likelihood does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.

#since experimental study, not testing the inclusion or exclusion of fixed effects, but rather assessing what interactions among fixed effects to include in final model

lmm.full <- lmer(days.forelimb ~ treatment + (1|clutchtank), data = devo.data, na.action = na.omit)
lmm.null <- lmer(days.forelimb ~ (1|clutchtank), data = devo.data, na.action = na.omit)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(lmm.full) #tests for over- and under-dispersion
testZeroInflation(lmm.full) #tests if more zeroes than expected
testCategorical(lmm.full, catPred = devo.data$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(lmm.full)

#based on DHARMaoutputs we have an issue with heteroskedasticity across treatment. A rule of thumb is that linear models are fairly robust to heterogeneity of variance so long as the maximum variance is no more than 4× greater than the minimum variance, which is demonstrated by the following ratio. So if this ratio is <4, you are likely okay with running the linear model
var(devo.data$days.forelimb[devo.data$treatment == "high density"], na.rm=T) / var(devo.data$days.forelimb[devo.data$treatment == "low density"], na.rm=T)

#but we can also check using a model that accounts for heteroskedasticity and allows the error variance to depend on fixed effects in addition to random effects to vary using nlme package; uses a generalized least squares regression (gls), which is a form of weighted regression
lme.full <- lme(fixed = days.forelimb ~ treatment,
                    random = (~1|clutchtank),
                    data = devo.data,
                    na.action = na.omit)

lme.full.var <- lme(fixed = days.forelimb ~ treatment,
                random = (~1|clutchtank),
                data = devo.data,
                weights = varIdent(form=~1|treatment),
                na.action = na.omit)

# check assumptions
par(mfrow = c(1,1))
DHARMa::testCategorical(lmm.full, catPred = devo.data$treatment)
plot(lme.full, resid(., type = "p") ~ fitted(.) | treatment, abline = 0 )
plot(lme.full.var, resid(., type = "p") ~ fitted(.) | treatment, abline = 0 )

# model selection using likelihood ratio test to determine if accounting for heterskedasticity significantly improves gls model fit
anova(lme.full, lme.full.var)

# Final Model: 
# check if different results if white.adjust set to true to allow for heteroskedasticity
car::Anova(lmm.full, white.adjust=TRUE, type = "II")
car::Anova(lmm.full, white.adjust=FALSE, type = "II")
summary(lmm.full)

final.mod = lme.full
summary(final.mod)


# ANALYZE DATA: Effect of rearing density and clutch:larval tank on larval duration for FIRST 10%---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but maximum likelihood does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.

#since experimental study, not testing the inclusion or exclusion of fixed effects, but rather assessing what interactions among fixed effects to include in final model

lmm.full <- lmer(days.forelimb ~ treatment + (1|clutchtank), data = devo.data.first10[devo.data.first10$first.10perc == "yes",], na.action = na.omit)
lmm.null <- lmer(days.forelimb ~ (1|clutchtank), data = devo.data.first10[devo.data.first10$first.10perc == "yes",], na.action = na.omit)

# check assumptions
simulateResiduals(lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(lmm.full) #tests for over- and under-dispersion
testZeroInflation(lmm.full) #tests if more zeroes than expected
testCategorical(lmm.full, catPred = devo.data.first10$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(lmm.full)

anova(lmm.full, lmm.null)

# Final Model: 
final.mod = lmm.full
car::Anova(final.mod, type = "II")
summary(final.mod)



# ANALYZE DATA: Effect of rearing density and clutch:larval tank on metamorphosis duration ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but maximum likelihood does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups.

#since experimental study, not testing the inclusion or exclusion of fixed effects, but rather assessing what interactions among fixed effects to include in final model
#for some reason need to subset out the dataset to not include NAs 
devo.data.temp = devo.data[is.na(devo.data$days.forelimb.tail) == FALSE,]
devo.data.temp = devo.data.temp[is.na(devo.data.temp$treatment) == FALSE,]

lmm.full <- lmer(days.forelimb.tail ~ treatment + (1|clutchtank), data = devo.data.temp, na.action = na.omit)

lmm.log <- lmer(log(days.forelimb.tail) ~ treatment + (1|clutchtank), data = devo.data.temp, na.action = na.omit)

lmm.pois <- glmer(days.forelimb.tail ~ treatment + (1|clutchtank), data = devo.data.temp, na.action = na.omit, family = poisson(link = "identity"))

lmm.pois.log <- glmer(days.forelimb.tail ~ treatment + (1|clutchtank), data = devo.data.temp, na.action = na.omit, family = poisson(link = "log"))

lmm.nb <- glmer(days.forelimb.tail ~ treatment + (1|clutchtank), data = devo.data.temp, na.action = na.omit, family=MASS::negative.binomial(theta=1.75))

lmm.null <- lmer(days.forelimb.tail ~ (1|clutchtank), data = devo.data.temp, na.action = na.omit)

# check assumptions
simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = lmm.log, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = lmm.pois, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = lmm.pois.log, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = lmm.nb, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.full, form = devo.data.temp$treatment)
plotResiduals(lmm.log, form = devo.data.temp$treatment)
plotResiduals(lmm.pois, form = devo.data.temp$treatment)
plotResiduals(lmm.pois.log, form = devo.data.temp$treatment)
plotResiduals(lmm.nb, form = devo.data.temp$treatment)

testDispersion(lmm.full)
testDispersion(lmm.pois) #tests for over- and under-dispersion
testCategorical(lmm.full, catPred = devo.data.temp$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.


#based on DHARMaoutputs we have an issue with heteroskedasticity across treatment. A rule of thumb is that linear models are fairly robust to heterogeneity of variance so long as the maximum variance is no more than 4× greater than the minimum variance, which is demonstrated by the following ratio. So if this ratio is <4, you are likely okay with running the linear model
var(devo.data.temp$days.forelimb.tail[devo.data.temp$treatment == "high density"], na.rm=T) / var(devo.data.temp$days.forelimb.tail[devo.data.temp$treatment == "low density"], na.rm=T)

#but we can also check using a model that accounts for heteroskedasticity and allows the error variance to depend on fixed effects in addition to random effects to vary using nlme package; uses a generalized least squares regression (gls), which is a form of weighted regression
lme.full <- lme(fixed = days.forelimb.tail ~ treatment,
                random = (~1|clutchtank),
                data = devo.data.temp,
                na.action = na.omit)

lme.full.var <- lme(fixed = days.forelimb.tail ~ treatment,
                    random = (~1|clutchtank),
                    data = devo.data.temp,
                    weights = varIdent(form=~1|treatment),
                    na.action = na.omit)

# model selection using likelihood ratio test to determine if accounting for heterskedasticity significantly improves gls model fit (it doesn't)
anova(lme.full, lme.full.var)

# Final Model: 
# check if different results if white.adjust set to true to allow for heteroskedasticity
car::Anova(lmm.full, white.adjust=TRUE, type = "II")
car::Anova(lmm.full, white.adjust=FALSE, type = "II")
summary(lmm.full)

anova(lmm.full, lmm.null)

final.mod = lmm.full
summary(final.mod)


# ANALYZE DATA: Effect of rearing density, mean days forelimb, and clutch on average gosner stage at end of experiment ---------------------
lmm.full <- lmer(tank.devo.stage ~ treatment*mean.days.forelimb + (1|clutch), data=survi.data.exp, na.action = na.omit)

lmm.nointxn <- lmer(tank.devo.stage ~ treatment + mean.days.forelimb + (1|clutch), data=survi.data.exp, na.action = na.omit)

lmm.null <- lmer(tank.devo.stage ~ (1|clutch), data=survi.data.exp, na.action = na.omit)
anova(lmm.full, lmm.nointxn, lmm.null)

# model selection using likelihood ratio test + checking whether random effect important to include
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
summary(lmm.nointxn)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn, quantreg=T, plot = T) #provides summary of model fitting tests

testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = survi.data.exp$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)


# ANALYZE DATA: Effect of rearing density, water level reduction, week, and larv tank id nested within clutch on % metamorphosing individuals weekly up until metamorphosis ---------------------

# model definition 
lmm.full <- lmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, weights = (seed.num-leth.samp.num.cumul))

glmm.full <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.probit <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial(link="probit"), weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.slopes <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc + (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn.slopes <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment + week + water.level.reduc + (week||clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))

glmm.null <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ (1|clutchtank), data=survi.data.tad.exp, na.action = na.omit, family = binomial, weights = (seed.num-leth.samp.num.cumul))


# check assumptions
simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
simulateResiduals(fittedModel = glmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
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

#examine the plots for for each tank, treatment, and clutch, regress metamorphosis proportion by week using linear and then using linear and you can see how much beter the glm is and how it does look like distinct slopes for high density vs low density
ggplot(survi.data.tad.exp, 
       aes(x = week, y = prop.seed.forelimb, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  #geom_smooth(data = survi.data.tad.exp, method = "lm", aes(x=week, y=prop.seed.forelimb), inherit.aes = F, se = F, color="black")
  geom_smooth(data = survi.data.tad.exp, method = "glm", method.args = list(family = "binomial"), aes(x=week, y=prop.seed.forelimb), inherit.aes = F, se = F, color="black")

ggplot() +
  stat_summary(data = survi.data.tad.exp, fun.y=mean, geom="line", size = 0.5, alpha = 0.5, aes(x = week, y = prop.seed.forelimb, color = treatment, group = clutchtank)) +
  stat_summary(data = survi.data.tad.exp, fun.y=mean, geom="point", color = "black", pch=21, alpha = 0.8, size=2, aes(x = week, y = prop.seed.forelimb, fill=treatment)) +
  #geom_smooth(data = survi.data.tad.exp, method = "lm", aes(x=week, y=prop.seed.forelimb), inherit.aes = F, se = F, color="black")
  geom_smooth(data = survi.data.tad.exp, method = "glm", method.args = list(family = "binomial"(link = "logit")), aes(x=week, y=prop.seed.forelimb, color = treatment, group = treatment, linetype = treatment), inherit.aes = F, se = F, size = 1.5) +
  scale_x_continuous(limits = c(1,11), breaks = seq(1,11,1))

# model selection using AIC among the logit link models
anova(glmm.full, glmm.full.slopes, glmm.nointxn, glmm.nointxn.slopes)

#but assumption plots indicate that we have a problem with overdispersion, which occurs when error (residuals) are more variable than expected from the theorized distribution
glmm.full.slopes.quasi <- MASS::glmmPQL(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ treatment*week + water.level.reduc,
              random=~week|clutchtank,
              family="quasibinomial",
              weights = (seed.num-leth.samp.num.cumul),
              data=survi.data.tad.exp)
#having trouble converging

#trying beta bin
survi.data.tad.exp$survi.num = survi.data.tad.exp$seed.num-survi.data.tad.exp$leth.samp.num.cumul
glmm.full.beta <- aod::betabin(formula = cbind(metamorph.num.cumul, survi.num - metamorph.num.cumul) ~ treatment*week + water.level.reduc, 
             random = ~clutchtank, 
             data = survi.data.tad.exp, 
             link = "logit", 
            warnings = FALSE, na.action = na.omit, fixpar = list(),
        hessian = TRUE, control = list(maxit = 2000))
#this also isn't working

#another option to try is modeling (0,1) data to see if that helps TRY MODELING A POLY2 OR POLY3 FUNCTION! BUT RENAME SURVI.DATA.TAD.EXP AS MORPH.DATA.TAD.EXP.LONGFORM, SINCE THAT'S MORE REALISTIC. Note some models estimate intercept and slope, separately, by random factor: (1 | random.factor) + (0 + fixed.factor | random.factor). An alternative way to write this is using the double-bar notation fixed.factor + (fixed.factor || random.factor). when graphing the data, while it looks like a week^2 or week^3 could be helpful, this is likely just because we lose clutches b and c for week 11 so the mean is artificially brought down...because of this, keeping the formula linear.
glmm.full.01 <- glmer(status ~ treatment*week + water.level.reduc + (1|clutchtank), data=morph.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes <- glmer(status ~ treatment*week + water.level.reduc + (week||clutchtank), data=morph.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ treatment + week + water.level.reduc + (1|clutchtank), data=morph.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.slopes <- glmer(status ~ treatment + week + water.level.reduc + (week||clutchtank), data=morph.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.01 <- glmer(status ~ (1|clutchtank), data=morph.data.tad.exp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

anova(glmm.full.01, glmm.full.01.slopes, glmm.nointxn.01, glmm.nointxn.01.poly2, glmm.nointxn.01.slopes, glmm.null.01)

#double-check assumptions

simulateResiduals(fittedModel = glmm.full.01.slopes, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01.slopes, form = morph.data.tad.exp.longform$treatment)
plotResiduals(glmm.full.01.slopes, form = morph.data.tad.exp.longform$week)
plotResiduals(glmm.full.01.slopes, form = morph.data.tad.exp.longform$water.level.reduc)
testDispersion(glmm.full.01.slopes)

# Final Model: although poly3 better supported, it causes a VERY weird looking graph so I think that actually the cubic root transformation is not an appropriate assumption for this
final.mod = glmm.full.01.slopes
summary(final.mod)

# create dataframe of predicted values that can be plotted on ggplot later
predicted.df.prop.morph <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, week, and larv tank id nested within clutch on % metamorphosing individuals CALCULATED ONLY FROM SURVIVING INDIVIDUALS weekly up until metamorphosis ---------------------

# model definition 
glmm.full.01<- glmer(status ~ treatment*week + water.level.reduc + (1|clutchtank), data=morph.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes<- glmer(status ~ treatment*week + water.level.reduc + (week||clutchtank), data=morph.data.tad.exp.longformsurvi, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tad.exp.longformsurvi$treatment)
plotResiduals(glmm.full.01, form = survi.data.tad.exp.longformsurvi$week)
plotResiduals(glmm.full.01, form = survi.data.tad.exp.longformsurvi$water.level.reduc)
testDispersion(glmm.full.01)
testOutliers(glmm.full.01)

anova(glmm.full.01, glmm.full.01.slopes)

#final model - allowing random slopes of treatment across week is not supported
summary(glmm.full.01)


# PLOT DATASETS: Effect of water level reduction on larval duration -----------------------

# plotted with clutch separated under each density group 
ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = water.level.reduc, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
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
  scale_x_discrete(name = "was water level reduced?")



# PLOT DATASETS: Effect of rearing density on larval duration -----------------------

# plotted with clutch separated under each density group 
plot.devo.1 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = clutch, color = treatment, show.legend = FALSE)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days)") +
  scale_x_discrete(name = "clutches separated")

# plotted with clutch lumped together under each density group for FIRST SIX
plot.devo.2 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = treatment, color = treatment, show.legend = FALSE)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days) - first six") +
  scale_x_discrete(name = "clutches combined")

plot.devo.2 <- ggMarginal(plot.devo.2, margins = "y", groupColour = TRUE, groupFill = TRUE)


# plotted with clutch lumped together under each density group for FIRST 10%
plot.devo.2first10perc <- ggplot(data = devo.data.first10[devo.data.first10$first.10perc == "yes",], aes(y=days.forelimb, x = treatment, color = treatment, show.legend = FALSE)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7, show.legend = FALSE) +
  geom_boxplot(alpha = 0.75, size = 0.75, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days) - first 10%") +
  scale_x_discrete(name = "clutches combined")

plot.devo.2first10perc <- ggMarginal(plot.devo.2first10perc, margins = "y", groupColour = TRUE, groupFill = TRUE)


# PLOT DATASETS: Effect of rearing density on proportion tank initiated metamorphosis by experiment completion ---------------------
# x-y plot with summarized mean and +/- 1 se for all metrics

# percent metamorphosed from tank seeding to tank close-out without clutch included
ggplot(data = survi.data.tad.exp, aes(y=prop.seed.forelimb*100, x = week, color = treatment)) + 
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

# percent metamorphosed from tank seeding to tank close-out without clutch included but now including predicted fit from model and plotting all tanks across time
plot.devo.3 <- ggplot() + 
  geom_jitter(width=0.4, size = 2, alpha = 0.7, data = survi.data.tad.exp, pch = 21, aes(y=prop.seed.forelimb, x = week, fill = treatment, color = treatment)) +
  geom_ribbon(data = predicted.df.prop.morph,
              alpha = 0.5,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group, fill = group), show.legend = FALSE) +
  geom_line(data = predicted.df.prop.morph,
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
  scale_y_continuous(name = "percent metamorphosed") + #, limits = c(0,0.35), breaks = c(0,.1,.2,.3), labels = c(0,10,20,30)) +
  scale_x_continuous(name = "age (weeks)", breaks = seq(0, 12, by = 1))

# raw number metamorphosed from tank seeding to tank close-out without clutch included
ggplot(data = survi.data.tad.exp, aes(y=metamorph.num.cumul, x = week, color = treatment)) + 
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
                legend.position = "none",
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

# number metamorphosed from tank seeding to tank close-out without clutch but now plotting all tanks across time
plot.devo.4 <- ggplot(data = survi.data.tad.exp, aes(y=metamorph.num.cumul, x = week, color = treatment)) + 
  geom_point(size = 2, alpha = 0.7, data = survi.data.tad.exp, pch = 21, aes(y=metamorph.num.cumul, x = week, fill = treatment)) +
  geom_point(size = 2, alpha = 1, data = survi.data.tad.exp, pch = 21, aes(y=metamorph.num.cumul, x = week, color = treatment)) +
  geom_line(size = 0.5, alpha = 0.7, data = survi.data.tad.exp, aes(y=metamorph.num.cumul, x = week, color = treatment, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=1, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment)) +
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
  scale_y_continuous(name = "number metamorphosed") +
  scale_x_continuous(name = "age (weeks)", breaks = seq(0, 12, by = 1))


# PLOT DATASETS: Panel plot of little-to-no effect of rearing density on developmental metrics ---------------------
ggarrange(plot.devo.3, plot.devo.4, plot.devo.1, plot.devo.2,
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))
# show.legend = FALSE in boxplots so that the common legend pulls from the other plots

#ggarrange doesn't show density plot from ggmarginal so trying this way
gridExtra::grid.arrange(plot.devo.3, plot.devo.4, plot.devo.1, plot.devo.2, 
                        ncol=2, nrow=2)



# PLOT DATASETS: Effect of rearing density on gosner stage for individuals from each tank (includes all metamorphosed individuals static at GS 42) -----------------------
ggplot(data = survi.data.exp, aes(y=tank.devo.stage, x = treatment, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  geom_boxplot(alpha = 0.6, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mean tank gosner stage at end of experiment", limits = c(26, 42), breaks = seq(26,42,2)) +
  scale_x_discrete()

#histogram representing distribution of tank means
ggplot() + 
  geom_histogram(data = survi.data.exp, na.rm = TRUE, binwidth = 0.5, color = "black", alpha = 0.7, position = "identity", aes(x = tank.devo.stage, fill = treatment)) + 
  facet_grid(rows = vars(treatment)) +
  geom_density(data = survi.data.exp, na.rm = TRUE, alpha = 0.5, aes(x = tank.devo.stage, color = treatment, fill = treatment)) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,5)) +
  scale_x_continuous(name = "gosner stage at end of experiment", limits = c(26, 42), breaks = seq(26,42,2))

#histogram representing distributions not plotted by tank mean
ggplot() + 
  geom_histogram(data = data.gstage, na.rm = TRUE, binwidth = 1, color = "black", alpha = 0.7, position = "identity", aes(x = devo.stage, fill = treatment)) + 
  facet_grid(rows = vars(treatment)) +
  geom_density(data = data.gstage, na.rm = TRUE, alpha = 0.5, aes(x = devo.stage, color = treatment, fill = treatment)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(limits = c(0,50)) +
  scale_x_continuous(name = "gosner stage at end of experiment", limits = c(26, 42), breaks = seq(26,42,2))


# PLOT DATASETS: Effect of rearing density, average gosner stage on size for individuals from each tank (includes all metamorphosed individuals static at GS 42) -----------------------
ggplot(data = survi.data.exp, aes(y=mean.mm.mass.g, x = tank.devo.stage, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mean mass at metamorphosis for first 6 individuals") +
  scale_x_continuous(name = "mean tank gosner stage at end of experiment", limits = c(26, 42), breaks = seq(26,42,2))


ggplot(data = survi.data.exp, aes(y=tank.devo.stage, x = mean.mm.mass.g, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mean tank gosner stage at end of experiment", limits = c(26, 42), breaks = seq(26,42,2)) +
  scale_x_continuous(name = "mean mass at metamorphosis", limits = c(0.1, 0.25))


# SUMMARY TABLES: Developmental metrics ---------------------

#reset working directory to be outputs folder
setwd("~/Desktop/R Working Directory/Outputs")

#tadpole summary
devo.data.mm.summ = devo.data %>%
  group_by(treatment) %>%
  filter(first.six == "yes") %>%
  summarise(mean.days.forelimb = round(mean(days.forelimb, na.rm = TRUE), 1),
            sd.days.forelimb = round(sd(days.forelimb, na.rm = TRUE), 1),
            min.days.forelimb = round(min(days.forelimb, na.rm = TRUE), 1),
            max.days.forelimb = round(max(days.forelimb, na.rm = TRUE), 1),
            range.days.forelimb = round(max(days.forelimb, na.rm = TRUE) - min(days.forelimb, na.rm = TRUE), 1),
            mean.days.completion = round(mean(days.forelimb.tail, na.rm = TRUE), 1),
            sd.days.completion = round(sd(days.forelimb.tail, na.rm = TRUE), 1)
  )

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


# SUMMARY TABLES: Developmental metrics with water level reduction included ---------------------

devo.data.summ = devo.data %>%
  group_by(treatment, water.level.reduc) %>%
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
  group_by(treatment, water.level.reduc) %>%
  summarise(num.seeded = sum(seed.num - leth.samp.num.cumul, na.rm = TRUE), 
            num.survi = sum(photo.num + metamorph.num.cumul, na.rm = TRUE),
            num.metamorph = sum(metamorph.num.cumul, na.rm = TRUE),
            mean.metamorph = round(mean(metamorph.num.cumul, na.rm = TRUE), 2),
            sd.metamorph = round(sd(metamorph.num.cumul, na.rm = TRUE), 2),
            mean.prop.seeded = round(mean(prop.seed.forelimb, na.rm = TRUE), 3),
            sd.prop.seeded = round(sd(prop.seed.forelimb, na.rm = TRUE), 3),
            mean.prop.survi = round(mean(prop.surv.forelimb, na.rm = TRUE), 3),
            sd.prop.survi = round(sd(prop.surv.forelimb, na.rm = TRUE), 3)
  )