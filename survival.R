# README -------------------------
# this script (1) compiles (cleans, joins) databases to be used in subsequent analyses, (2) creates figures of the survivorship results for the 2025 ranid developmental paper

# LOAD PACKAGES and SET WD ------------------------------------
rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(coxme) #don't need this package if not running Cox proportional hazards models
library(car)
library(lme4)
library(DHARMa)
library(MuMIn)
library(RColorBrewer)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for survival data and development data
survi.data.tad = read.csv("Database_Survivorship - Tadpole Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.exp = read.csv("Database_Survivorship - Tadpole Experiment Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
survi.data.juv = read.csv("Database_Survivorship - Froglet_Toadlet Weekly Survivorship.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database

morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
survi.data.tad = survi.data.tad[survi.data.tad$gs.code == "RS",]
survi.data.exp = survi.data.exp[survi.data.exp$gs.code == "RS",]
survi.data.juv = survi.data.juv[survi.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS" & is.na(devo.data$gs.code) == FALSE,]
morph.data.mm = morph.data.mm[morph.data.mm$gs.code == "RS",]

# subset databases to only include weeks up to 18
survi.data.juv = survi.data.juv[is.na(survi.data.juv$postmm.week) == FALSE,]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
survi.data.tad$treatment[survi.data.tad$treatment == "control"] = "low density"
survi.data.exp$treatment[survi.data.exp$treatment == "control"] = "low density"
survi.data.juv$treatment[survi.data.juv$treatment == "control"] = "low density"
devo.data$treatment[devo.data$treatment == "control"] = "low density"
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
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")

#remove any NAs that have occurred
devo.data = devo.data[is.na(devo.data$gs) == FALSE,]

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


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.mm -----------------------
#create unique.id.indiv
morph.data.mm$unique.id.indiv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, morph.data.mm$animal.id.num, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$unique.id, devo.data$animal.id, sep = "_")

# create column to store the developmental data
morph.data.mm$days.forelimb = NA

# fill developmental data column for morph.data.mm
for(i in 1:length(unique(morph.data.mm$unique.id.indiv))){
  morph.data.mm$days.forelimb[morph.data.mm$unique.id.indiv == unique(morph.data.mm$unique.id.indiv)[i]] <- devo.data$days.forelimb[devo.data$unique.id.indiv ==  unique(morph.data.mm$unique.id.indiv)[i]]
}


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
#jitter points
survi.data.tad.exp$prop.survi.j = jitter(survi.data.tad.exp$prop.survi, amount = 0.001)

fig.7a <- ggplot() + 
  
  #individual tanks
  geom_point(size = 1, alpha = 0.8, data = survi.data.tad.exp, pch = 21, aes(y=prop.survi.j , x = as.numeric(week), color = treatment, fill = treatment)) +
  geom_line(size = 0.5, alpha = 0.7, data = survi.data.tad.exp, aes(y=prop.survi.j , x = as.numeric(week), color = treatment, group = clutchtank)) +
  
  #treatment means
  stat_summary(data = survi.data.tad.exp, fun=mean, geom="line", size = 0.8, color = "black", aes(x = as.numeric(week), y = prop.survi, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = survi.data.tad.exp, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = as.numeric(week), y = prop.survi, group = treatment)) +
  stat_summary(data = survi.data.tad.exp, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = as.numeric(week), y = prop.survi, fill=treatment), show.legend = TRUE) + 
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = -11.29", sep = ""))), # sum of coefficient for intercept, high density, week, and high density*week
           x=1, y = 0.88, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 11.59", sep = ""))),
           x=1, y = 0.87, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.001",
           x=1, y = 0.86, color = "black", size = 5, hjust=0) +
  
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

  scale_y_continuous(name = "percent larval survivorship", limits=c(0.85,1.001), breaks = seq(0.85,1.0,0.05), labels = seq(85,100,5)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


# PLOT DATASETS: Effect of predictors on survival after metamorphosis ---------------------

# alternative showing weird interaction with low density
temp2 = ggpredict(final.mod, terms = c("postmm.week [1,9,18]", "treatment", "mean.mm.mass.g [all]", "mean.days.forelimb [all]"), type = "random", interval = "confidence")
colnames(temp2)[6:8] = c("treatment", "mean.mm.mass.g", "mean.days.forelimb")
colnames(temp2)[1] = "postmm.week"
temp2$treatment = factor(temp2$treatment)
temp2$mean.mm.mass.g = as.character(temp2$mean.mm.mass.g)
temp2$mean.mm.mass.g = as.numeric(temp2$mean.mm.mass.g)

fig.6c <- ggplot() + 
  facet_grid(cols=vars(as.factor(postmm.week)), rows = vars(treatment)) + 
  
  geom_point(alpha = 0.7, pch = 21, color = "black",
             data = survi.data.juv[survi.data.juv$postmm.week == 1 | survi.data.juv$postmm.week == 9 | survi.data.juv$postmm.week == 18,], 
             aes(y=prop.survi, x = mean.mm.mass.g, fill = mean.days.forelimb, size = mean.mm.mass.g), show.legend = TRUE) +

  #predicted from model: at longer larval duration, the negative effect of metamorphic mass on larval survivorship reverses, such that smaller, longer larval duration tanks show better survivorship BUT ONLY FOR HIGH DENSITY 
  geom_line(data=temp2,
            aes(x = mean.mm.mass.g, y=predicted, color=mean.days.forelimb, group=mean.days.forelimb),
            show.legend = FALSE) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_fill_gradient(low = "gray85", high = "gray1") +
  scale_size(range = c(1.5, 6), limits = c(0.132,0.192), breaks = seq(0.132,0.192, 0.02)) + #, limits = c(0.13,0.191), breaks = seq(0.13, 0.19, 0.02), labels = c("0.13", "0.15", "0.17", "0.19")) +
  labs(fill = "tank-level mean larval duration (days)", size = "tank-level mean metamorphic mass (g)") +
  
  theme_bw() +
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=20, angle = 45, hjust= 1, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text( #increases size of labels in facet grid
          size = 20)) +
  scale_y_continuous(name = "percent juvenile survivorship", labels = c(0,25,50,75,100)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "tank-level mean metamorphic mass (g)", limits = c(0.13,0.2), breaks = seq(0.13,0.2,0.02), labels = c("0.13", "0.15", "0.17", "0.19"))

# add low density and high density colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
g.c <- ggplot_gtable(ggplot_build(fig.5c))
strip_r <- which(grepl('strip-r', g.5c$layout$name))
fills <- c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])
colors <- c("grey10", "white")
k <- 1
for(i in strip_r){
  j <- which(grepl('rect', g.5c$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.5c$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  g.5c$grobs[[i]]$grobs[[1]]$children[[2]]$children[[1]]$gp$col <- colors[k] #only colors the row labels white
  k <- k+1
}
grid::grid.draw(g.5c) #produces final plot with colored headers only for rows

  
#create dataframe with rectangle dimensions to block annotation on one facet but not the other
rect<-data.frame(treatment = c("low density", "high density"), 
                 xmin = c(20, 22), xmax = c(20.5, 31), 
                 ymin = c(0.6, 0.8), ymax = c(0.65, 1.0), 
                 alpha = c(1, 1),
                 fill = c("green", "green"))
rect$treatment = factor(rect$treatment)

fig.7b <- ggplot() + 
  facet_grid(cols=vars(treatment)) + 

  stat_summary(data = survi.data.juv, fun = mean,
               geom = "errorbar", position = position_jitter(width = 0.35, height = 0.0025, seed = 0.1),
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, alpha=0.7, aes(x = postmm.week, y = prop.survi, group = clutchtank), show.legend = F) +
  
  stat_summary(data = survi.data.juv, fun = mean,
               geom = "line", size = 0.8, position = position_jitter(width = 0.35, height = 0.0025, seed = 0.1),
               alpha=0.75, aes(x = postmm.week, y = prop.survi, group = clutchtank, color = mean.mm.mass.g), show.legend = F) +
  
  stat_summary(data = survi.data.juv, fun = mean,
               geom = "point", pch = 21, position = position_jitter(width = 0.35, height = 0.0025, seed = 0.1), color = "black", stroke = 0.5, size = 3,
               alpha=0.75, aes(x = postmm.week, y = prop.survi, group = clutchtank, fill = mean.mm.mass.g), show.legend = T) +
  
  scale_color_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  scale_fill_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = -1.95", sep = ""))), # sum of coefficient for intercept, mean larval duration, mean metamorphic mass, post-metamorphic age, and high density*larval duration*post-metamorphic age
           x=24, y = 0.98, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 27.92", sep = ""))),
           x=24, y = 0.90, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.001",
           x=24, y = 0.82, color = "black", size = 5, hjust=0) +
  
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white",
            data = rect,
            inherit.aes = FALSE) +
  
  labs(fill = "tank mean metamorphic mass (g)") +
  
  theme_bw() +
  theme(legend.position = c(0.98,0.98),
        legend.justification = c("right", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text( #increases size of labels in facet grid
          size = 20, color = "white")) +
  scale_y_continuous(name = "percent juvenile survivorship", labels = c(0,25,50,75,100), breaks = seq(0,1,0.25), limits = c(-0.008,1.0025)) + #binomial can only be fit 0-1 but relabeled to be percent rather than proportion survival
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = seq(1, 31, by = 5), expand=c(0.02, 0), limits = c(0.5, 31)) #expand slightly widens x-axis


# add low density and high density colors to facet grid (code from: https://github.com/tidyverse/ggplot2/issues/2096)
g.7b <- ggplot_gtable(ggplot_build(fig.7b))
strip_t <- which(grepl('strip-t', g.7b$layout$name))
fills <- c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])
k <- 1
for(i in strip_t){
  j <- which(grepl('rect', g.7b$grobs[[i]]$grobs[[1]]$childrenOrder))
  g.7b$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid::grid.draw(g.7b) #produces final plot with colored headers




# PLOT DATA: Create panel plot with all survival data across all sampling points --------------
ggarrange(fig.a, 
          ggarrange(fig.3b, fig.3c,
                    ncol = 2,
                    nrow = 1,
                    labels = c("b", "c"),
                    common.legend = FALSE,
                    legend = "bottom",
          font.label = list(size = 20, color = "black")),
          ggarrange(fig.3d,
                    ncol=1,
                    nrow=1,
                    labels = c("d"),
                    legend = NULL),
          ncol = 1,
          nrow = 3,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a"),
          font.label = list(size = 20, color = "black"))

cowplot::plot_grid(grid::grid.draw(g.5b),
          grid::grid.draw(g.5c),
          ncol = 1,
          nrow = 3)

## I think I need to set the heights of the grid.draw figures manually before feeding into grid.arrange
gridExtra::grid.arrange(fig.5a,
                        gridExtra::grid.arrange(grid::grid.draw(g.5c), grid::grid.draw(g.5b),
                                     ncol = 1, nrow = 2),
                        ncol=2, nrow=1,
                        )

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

# PANEL PLOTS + EXPORT PLOTS: Create panel plot with mass as a function of treatment and larval duration across all sampling points --------------

png("~/Desktop/R Working Directory/Plots/Figure7.png", units = "in", res = 300, width = 16, height = 12)
ggarrange(fig.7a, g.7b, 
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()


# DOWNLOAD DATA FOR JOURNAL SUBMISSION  ---------------------
#set wd
setwd("~/Desktop/R Working Directory/2023_ranid_devo_plasticity/Submission_RoySocB")

write.csv(survi.data.tad[,c(3:11,16:21, 25:26)], file = "Database_Survivorship_Tadpole_Weekly.csv", row.names = FALSE)

write.csv(survi.data.juv[,c(5, 7, 10:14, 16:17, 19:21, 23:25, 27, 30)], file = "Database_Survivorship_Juvenile_Weekly.csv", row.names = FALSE)

write.csv(survi.data.juv.longform, file = "Database_Survivorship_Juvenile_Longform.csv", row.names = FALSE)

