# README -------------------------
# this script (1) compiles (cleans, joins) databases to be used in subsequent analyses, (2) creates figures of the developmental results for the 2025 ranid developmental paper

# LOAD PACKAGES and SET WD ------------------------------------
rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)
library(DHARMa)
library(ggsignif)
library(MuMIn)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for developmental data and survival data (provides tank-level percent metamorphosed by end of experiment)
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data.nonmm = read.csv("Database_Metamorphosis - Non-Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database

survi.data.exp = read.csv("Database_Survivorship - Tadpole Experiment Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.tad = read.csv("Database_Survivorship - Tadpole Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
devo.data = devo.data[devo.data$gs.code == "RS",]
devo.data.nonmm = devo.data.nonmm[devo.data.nonmm$gs.code == "RS",]
survi.data.exp = survi.data.exp[survi.data.exp$gs.code == "RS",]
survi.data.tad = survi.data.tad[survi.data.tad$gs.code == "RS",]

# subset databases to only include first six individuals and non-overflow
devo.data = devo.data[devo.data$treatment != "overflow",]
devo.data = devo.data[devo.data$first.six == "yes",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data$treatment[devo.data$treatment == "control"] = "low density"
devo.data.nonmm$treatment[devo.data.nonmm$treatment == "control"] = "low density"
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

devo.data.nonmm$treatment = factor(devo.data.nonmm$treatment, levels = c("low density", "high density"))
devo.data.nonmm$larv.tank.id = factor(devo.data.nonmm$larv.tank.id)

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
devo.data.nonmm$unique.id = paste(devo.data.nonmm$gs.code, devo.data.nonmm$clutch, devo.data.nonmm$larv.tank.id, sep = "_")
survi.data.tad$unique.id = paste(survi.data.tad$gs.code, survi.data.tad$clutch, survi.data.tad$larv.tank.id, survi.data.tad$week, sep = "_")

# rename two columns to match survi.data.tad for later binding
colnames(survi.data.exp)[11] = "photo.num"
colnames(survi.data.exp)[12] = "leth.samp.num.cumul"

#create column to calculate average gosner stage
devo.data.nonmm$avg.gosner.stage = (devo.data.nonmm$min.gosner.stage + devo.data.nonmm$max.gosner.stage)/2


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

survi.data.tad.exp$unique.id.week = paste(survi.data.tad.exp$unique.id, survi.data.tad.exp$week, sep = "_")


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
  survi.data.exp$tank.devo.stage[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] = mean(temp, na.rm = TRUE)
}


# COMPILE DATASETS: Create long.form Gosner staging INCLUDING METAMORPHOSED INDIVIDUALS  -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank

# fill developmental data column for survi.data.exp.
for(i in 1:length(unique(survi.data.exp$unique.id))){
  
  temp1 = devo.data.nonmm[devo.data.nonmm$unique.id == unique(survi.data.exp$unique.id)[i], c(1:6,13)]
  
  #create vector that includes all the ones that already metamorphosed as stage 42 + non-metamorphosed stage
  temp.mm = survi.data.exp$metamorph.num.cumul[survi.data.exp$unique.id == unique(survi.data.exp$unique.id)[i]] #number of individuals "stuck" at 42
  
  temp1[nrow(temp1) + temp.mm,] <- NA #add number of empty rows needed

  temp1[is.na(temp1[,1]) == TRUE,1] = temp1[1,1]
  temp1[is.na(temp1[,2]) == TRUE,2] = temp1[1,2]
  temp1[is.na(temp1[,3]) == TRUE,3] = temp1[1,3]
  temp1[is.na(temp1[,4]) == TRUE,4] = temp1[1,4]
  temp1[is.na(temp1[,5]) == TRUE,5] = temp1[1,5]
  temp1[is.na(temp1[,6]) == TRUE,6] = temp1[1,6]

  temp1[is.na(temp1[,7]) == TRUE,7] = 42
  
  if(i==1){devo.data.longform = temp1}else{
    devo.data.longform = rbind(devo.data.longform, temp1)
  }
  
rm(temp1, temp.mm)

}

devo.data.longform$clutchtank = factor(paste(devo.data.longform$clutch, devo.data.longform$larv.tank.id, sep = "_"))


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
morph.data.tad.exp.longform$clutchtank = factor(paste(morph.data.tad.exp.longform$clutch, morph.data.tad.exp.longform$larv.tank.id, sep = "_"))
morph.data.tad.exp.longformsurvi$clutchtank = factor(paste(morph.data.tad.exp.longformsurvi$clutch, morph.data.tad.exp.longformsurvi$larv.tank.id, sep = "_"))
survi.data.tad.exp$clutchtank = factor(paste(survi.data.tad.exp$clutch, survi.data.tad.exp$larv.tank.id, sep = "_"))
survi.data.exp$clutchtank = factor(paste(survi.data.exp$clutch, survi.data.exp$larv.tank.id, sep = "_"))

data.gstage$clutchtank = factor(paste(data.gstage$clutch, data.gstage$larv.tank.id, sep = "_"))



# COMPILE DATASETS: Create 10% dataset to identify the first 10% to metamorphose for each rearing density -----------------------
devo.data.first10 = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
devo.data.first10 = filter(devo.data.first10, gs.code == "RS")

# subset databases to only include non-overflow
devo.data.first10 = filter(devo.data.first10, treatment != "overflow")

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data.first10$treatment[devo.data.first10$treatment == "control"] = "low density"

#subset out only first 10% (first 4 for low density, first 10 for high density)
devo.data.first10 = filter(devo.data.first10, first.10perc == "yes")

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

#how many high density tanks did we remove that hadn't reached 10?
first10.summ = devo.data.first10 %>%
  group_by(treatment, clutch, larv.tank.id) %>%
  summarize(n = n())
#nrow(temp[first10.summ$treatment == "high density" & first10.summ$n < 10,]) #this provides number of tanks that we are removing by opportunistically looking at first 10% of tanks rather than first six individuals

#create clutchtank for later analysis
devo.data.first10$clutchtank = paste(devo.data.first10$clutch, devo.data.first10$larv.tank.id, sep="_")

# PLOT DATASETS: Effect of rearing density on larval duration -----------------------

# plotted with clutch lumped together under each density group for FIRST SIX
ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = treatment, color = treatment, show.legend = FALSE)) + 
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
fig.1c <- ggplot(data = filter(devo.data, first.10perc == "yes"), aes(y=days.forelimb, x = treatment, color = treatment, show.legend = FALSE)) + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), size = 1.5, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(alpha = 0.6, size = 1, show.legend = FALSE, aes(fill = treatment), outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  geom_signif(comparisons = list(c("low density", "high density")), 
              map_signif_level=TRUE,
              annotations = "p < 0.001",
              size = 1,
              tip_length = 0,
              vjust = -0.3,
              textsize = 6,
              color = "black") + 
  
  scale_y_continuous(name = "larval duration (days)", limits = c(43.75,76), breaks = seq(45,76,5), labels = seq(45,76,5)) +
  scale_x_discrete(name = "")


# PLOT DATASETS: Effect of rearing density on proportion tank initiated metamorphosis by experiment completion ---------------------
# x-y plot with summarized mean and +/- 1 se for all metrics

# percent metamorphosed from tank seeding to tank close-out without clutch included but now including predicted fit from model and plotting all tanks across time
fig.1a <- ggplot() + 
  
  #individual tanks
  #geom_jitter(width=0.4, size = 2, alpha = 0.7, data = survi.data.tad.exp, pch = 21, aes(y=prop.seed.forelimb, x = as.numeric(week), fill = treatment, color = treatment)) +
  geom_point(size = 1.5, alpha = 0.7, data = survi.data.tadexp, pch = 21, aes(y=prop.seed.forelimb, x = as.numeric(week), color = treatment, fill =treatment)) +
  geom_line(size = 0.5, alpha = 0.7, data = survi.data.tadexp, aes(y=prop.seed.forelimb, x = as.numeric(week), color = treatment, group = clutchtank)) +
  
  #treatment means
  stat_summary(data = survi.data.tadexp, fun=mean, geom="line", size = 1, color = "black", aes(x = as.numeric(week), y = prop.seed.forelimb, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = survi.data.tadexp, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = as.numeric(week), y = prop.seed.forelimb, group = treatment)) +
  stat_summary(data = survi.data.tadexp, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = as.numeric(week), y = prop.seed.forelimb, fill=treatment), show.legend = TRUE) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "% metamorphosed", limits = c(0,0.3), breaks = seq(0,0.3,0.1), labels = c(0,10,20,30)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 11, by = 1))


# raw number metamorphosed from tank seeding to tank close-out without clutch included and plotting all tanks across time with means overlaid
fig.1b <- ggplot(data = survi.data.tadexp, aes(y=metamorph.num.cumul, x = as.numeric(week))) + 
  
  #individual tanks
  #geom_jitter(width=0.4, size = 2, alpha = 0.7, data = survi.data.tad.exp, pch = 21, aes(y=prop.seed.forelimb, x = as.numeric(week), fill = treatment, color = treatment)) +
  geom_point(size = 1.5, alpha = 0.7, pch = 21, aes(color = treatment, fill = treatment)) +
  geom_line(size = 0.5, alpha = 0.7, aes(color = treatment, group = clutchtank)) +
  
  #treatment means
  stat_summary(fun=mean, geom="line", size = 1, color = "black", aes(group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(group = treatment)) +
  stat_summary(fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(fill=treatment), show.legend = TRUE) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "# metamorphosed", limits = c(0,20)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(1, 11, by = 1))


# PLOT DATASETS: Effect of rearing density on gosner stage for individuals from each tank (includes all metamorphosed individuals static at GS 42) -----------------------
fig.1d <- ggplot(data = survi.data.exp, aes(y=tank.devo.stage, x = treatment, color = treatment), show.legend = FALSE) + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), size = 1.5, alpha = 0.7, show.legend = FALSE, aes(group=clutchtank)) +
  
  geom_boxplot(alpha = 0.6, size = 1, show.legend = FALSE, aes(fill = treatment), outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  geom_signif(comparisons = list(c("low density", "high density")), 
              map_signif_level=TRUE,
              annotations = "p < 0.001",
              size = 1,
              tip_length = 0,
              vjust = -0.3,
              textsize = 6,
              color = "black") + 
  
  scale_y_continuous(name = "mean Gosner stage", limits = c(31.75,39.4), breaks = seq(32,39,1), labels = seq(32,39,1)) +
  scale_x_discrete(name = "")


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


# PANEL PLOTS + EXPORT PLOTS: Panel plot of developmental results ---------------------
pdf("~/Desktop/R Working Directory/Plots/Figure1.pdf", width = 11, height = 8.5)
ggarrange(fig.1a, fig.1b, fig.1c, fig.1d,
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"),
          vjust = 1.5)
dev.off()

png("~/Desktop/R Working Directory/Plots/Figure1.png", units = "in", res = 300, width = 8.5, height = 8.5)
ggarrange(fig.1a, fig.1b, fig.1c, fig.1d,
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"),
          vjust = 1.5)
dev.off()
# show.legend = FALSE in boxplots so that the common legend pulls from the other plots

#ggarrange doesn't show density plot from ggmarginal so trying this way
gridExtra::grid.arrange(fig.1a, fig.1b, 
                        ggExtra::ggMarginal(fig.1c, margins = "y", groupColour = TRUE, groupFill = TRUE), 
                        ggExtra::ggMarginal(fig.1d, margins = "y", groupColour = TRUE, groupFill = TRUE), 
                        ncol=2, nrow=2,
                        
                        )


# DOWNLOAD DATA FOR JOURNAL SUBMISSION  ---------------------
#set wd
setwd("~/Desktop/R Working Directory/2023_ranid_devo_plasticity/Submission_RoySocB")

write.csv(devo.data.first10[devo.data.first10$first.10perc == "yes", c(1:8, 9:10, 13:16, 30, 32)], file = "Database_Development.csv", row.names = FALSE)

write.csv(survi.data.exp[,c(1:9, 11:16, 18:19, 21:22)], file = "Database_Survivorship_Tadpole_Experiment.csv", row.names = FALSE) #need to do this here because includes gosner stage data in addition to survivorship data

write.csv(survi.data.tad.exp[,c(1:16,18)], file = "Database_Survivorship_Tadpole_WeeklyAndExperiment.csv", row.names = FALSE)

write.csv(morph.data.tad.exp.longform, file = "Database_Survivorship_Tadpole_WeeklyAndExperiment_Longform.csv", row.names = FALSE)

write.csv(morph.data.tad.exp.longformsurvi, file = "Database_Survivorship_Tadpole_WeeklyAndExperiment_LongformSurvi.csv", row.names = FALSE)


# SUMMARY TABLES: Developmental metrics ---------------------

#reset working directory to be outputs folder
setwd("~/Desktop/R Working Directory/Outputs")

#by week summary
devo.data.mm.summ.byweek = survi.data.tad.exp %>%
  group_by(treatment, week) %>%
  filter(gs.code== "RS") %>%
  summarise(mean.prop.forelimb = round(mean(prop.seed.forelimb, na.rm = TRUE), 2),
            sd.prop.forelimb = round(sd(prop.seed.forelimb, na.rm = TRUE), 2),
            min.prop.forelimb = round(min(prop.seed.forelimb, na.rm = TRUE), 2),
            max.prop.forelimb = round(max(prop.seed.forelimb, na.rm = TRUE), 2),
            range.prop.forelimb = round(max(prop.seed.forelimb, na.rm = TRUE) - min(prop.seed.forelimb, na.rm = TRUE), 2),
  )


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
  summarise(num.seeded = sum(seed.num - leth.samp.num.cumul, na.rm = TRUE), 
            num.survi = sum(photo.num + metamorph.num.cumul, na.rm = TRUE),
            num.metamorph = sum(metamorph.num.cumul, na.rm = TRUE),
            mean.metamorph = round(mean(metamorph.num.cumul, na.rm = TRUE), 2),
            sd.metamorph = round(sd(metamorph.num.cumul, na.rm = TRUE), 2),
            mean.prop.seeded = round(mean(prop.seed.forelimb, na.rm = TRUE), 3),
            sd.prop.seeded = round(sd(prop.seed.forelimb, na.rm = TRUE), 3),
            mean.prop.survi = round(mean(prop.surv.forelimb, na.rm = TRUE), 3),
            sd.prop.survi = round(sd(prop.surv.forelimb, na.rm = TRUE), 3),
            mean.tank.devo.stage = round(mean(tank.devo.stage, na.rm = TRUE), 1),
            sd.tank.devo.stage = round(sd(tank.devo.stage, na.rm = TRUE), 1)
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
  group_by(treatment) %>%
  summarise(num.seeded = sum(seed.num - leth.samp.num.cumul, na.rm = TRUE), 
            num.survi = sum(photo.num + metamorph.num.cumul, na.rm = TRUE),
            num.metamorph = sum(metamorph.num.cumul, na.rm = TRUE),
            min.prop.seeded = min(prop.seed.forelimb, na.rm = TRUE),
            max.prop.seed = max(prop.seed.forelimb, na.rm = TRUE),
            mean.metamorph = round(mean(metamorph.num.cumul, na.rm = TRUE), 2),
            sd.metamorph = round(sd(metamorph.num.cumul, na.rm = TRUE), 2),
            min.metamorph = min(metamorph.num.cumul, na.rm = TRUE),
            max.metamorph = max(metamorph.num.cumul, na.rm = TRUE),
            mean.prop.seeded = round(mean(prop.seed.forelimb, na.rm = TRUE), 3),
            sd.prop.seeded = round(sd(prop.seed.forelimb, na.rm = TRUE), 3),
            mean.prop.survi = round(mean(prop.surv.forelimb, na.rm = TRUE), 3),
            sd.prop.survi = round(sd(prop.surv.forelimb, na.rm = TRUE), 3)
  )