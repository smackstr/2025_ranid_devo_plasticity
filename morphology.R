# README -------------------------
# this script (1) compiles (cleans, joins) databases to be used in subsequent analyses, (2) creates figures of the morphological results for the 2025 ranid developmental paper

# LOAD PACKAGES and SET WD ------------------------------------

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)
library(lme4)
library(car)
library(DHARMa)
library(devtools)
library(stringr)
library(emmeans)
library(ggeffects)
library(MuMIn)
library(ggnewscale)
library(pwr)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for metamorph morphometrics, juvenile morphometrics, and developmental timing data
morph.data.tad = read.csv("Database_Morphometrics - Tadpole Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv = read.csv("Database_Morphometrics - Froglet_Toadlet Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data.nonmm = read.csv("Database_Metamorphosis - Non-Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") #read in metamorphosis timing log database

# subset databases to only include Rana sylvatica and remove NAs
morph.data.tad = morph.data.tad[morph.data.tad$gs.code == "RS",]
morph.data.mm = morph.data.mm[morph.data.mm$gs.code == "RS",]
morph.data.juv = morph.data.juv[morph.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS" & is.na(devo.data$gs.code) == FALSE,]

# subset tadpole database to only include hormones and morphometrics measurements
morph.data.tad = morph.data.tad[morph.data.tad$data.type == "hormones" | morph.data.tad$data.type == "morphometrics",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
morph.data.tad$treatment[morph.data.tad$treatment == "control"] = "low density"
morph.data.mm$treatment[morph.data.mm$treatment == "control"] = "low density"
morph.data.juv$treatment[morph.data.juv$treatment == "control"] = "low density"
devo.data$treatment[devo.data$treatment == "control"] = "low density"
devo.data.nonmm$treatment[devo.data.nonmm$treatment == "control"] = "low density"

# subset devo dataset to only include first six individuals and non-overflow individuals
#devo.data = devo.data[devo.data$first.six == "yes",]
#devo.data = devo.data[devo.data$treatment != "overflow",]

# subset morph.data.mm & morph.data.juv to only include first six individuals and non-overflow individuals
morph.data.mm = morph.data.mm[morph.data.mm$first.six == "yes",]
morph.data.mm = morph.data.mm[morph.data.mm$treatment != "overflow",]
morph.data.juv = morph.data.juv[morph.data.juv$treatment != "overflow",]

# subset morph.data.tad to only include non-overflow individuals and non-extra individuals
morph.data.tad = morph.data.tad[morph.data.tad$treatment == "high density" | morph.data.tad$treatment == "low density" ,]
morph.data.tad = morph.data.tad[which(str_detect(morph.data.tad$larv.tank.id, "extra") == FALSE),]

#change column classes
morph.data.mm$first.six = factor(morph.data.mm$first.six, levels = c("yes", "no"))
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))

morph.data.tad$treatment = factor(morph.data.tad$treatment, levels = c("low density", "high density"))
morph.data.mm$treatment = factor(morph.data.mm$treatment, levels = c("low density", "high density"))
morph.data.juv$treatment = factor(morph.data.juv$treatment, levels = c("low density", "high density"))
devo.data$treatment = factor(devo.data$treatment, levels = c("low density", "high density"))

morph.data.tad$larv.tank.id = factor(morph.data.tad$larv.tank.id)
morph.data.mm$larv.tank.id = factor(morph.data.mm$larv.tank.id)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)

morph.data.mm$juv.tank.id = factor(morph.data.mm$juv.tank.id)
morph.data.juv$juv.tank.id = factor(morph.data.juv$juv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)

morph.data.mm$post.mm.weeks = factor(morph.data.mm$post.mm.weeks, ordered = TRUE, levels = c("0")) #set weeks as ordered factor
morph.data.juv$post.mm.weeks = factor(morph.data.juv$post.mm.weeks, ordered = TRUE, levels = c("1-2", "4-6", "5-7", "8-10", "11-12", "12-14", "16-18", "30-31")) #set weeks as ordered factor

morph.data.tad$tmw.cm = as.numeric(morph.data.tad$tmw.cm)

devo.data.nonmm$treatment = factor(devo.data.nonmm$treatment, levels = c("low density", "high density"))
devo.data.nonmm$larv.tank.id = factor(devo.data.nonmm$larv.tank.id)

# only keep certain weeks (0, 1-2, 4-6, 8-10, 12-14) and when individual not measured twice within that week. some weeks had an individual measured twice - for example, if the individual was measured for both morphometrics and ephys, but only one of these should be considered since it's a repeat measure on an (unknown) individual
morph.data.juv = morph.data.juv[morph.data.juv$post.mm.weeks != "5-7" & morph.data.juv$post.mm.weeks != "11-12",]
morph.data.juv = morph.data.juv[morph.data.juv$repeat.measure == "no",]

# remove rows for NA tanks and rename extra tanks
morph.data.juv = morph.data.juv[is.na(morph.data.juv$juv.tank.id) == FALSE,]
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J002extra"] = "J002"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J003extra"] = "J003"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J005extra"] = "J005"

# create unique tank id for all datasets
morph.data.tad$unique.id = paste(morph.data.tad$gs.code, morph.data.tad$clutch, morph.data.tad$larv.tank.id, sep = "_")
morph.data.mm$unique.id.larv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, sep = "_")
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")
morph.data.juv$unique.id = paste(morph.data.juv$gs.code, morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.juv = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")
devo.data.nonmm$unique.id = paste(devo.data.nonmm$gs.code, devo.data.nonmm$clutch, devo.data.nonmm$larv.tank.id, sep = "_")

# create unique individual id for devo.data and morph.data.mm, since these are the two we have individual-level data for
morph.data.mm$unique.id.indiv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, morph.data.mm$animal.id, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, devo.data$animal.id, sep = "_")

#remove any NAs that have occurred
devo.data = devo.data[is.na(devo.data$gs) == FALSE,]
morph.data.tad = morph.data.tad[is.na(morph.data.tad$gs) == FALSE,]

#create column to calculate average gosner stage
devo.data.nonmm$avg.gosner.stage = (devo.data.nonmm$min.gosner.stage + devo.data.nonmm$max.gosner.stage)/2

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
                                             "unique.id")], 
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
                                             "unique.id")])

# create column so can plot post.mm.weeks on numeric scale
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "0"] = 0
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "1-2"] = mean(c(1,2))
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "4-6"] = mean(c(4,6))
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "8-10"] = mean(c(8,10))
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "12-14"] = mean(c(12,14))
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "16-18"] = mean(c(16,18))
morph.data.mm.juv$post.mm.weeks.num[morph.data.mm.juv$post.mm.weeks == "30-31"] = mean(c(30,31))

# create column so can plot post.mm.weeks on numeric scale
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "1-2"] = mean(c(1,2))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "4-6"] = mean(c(4,6))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "8-10"] = mean(c(8,10))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "12-14"] = mean(c(12,14))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "16-18"] = mean(c(16,18))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "30-31"] = mean(c(30,31))


# COMPILE DATASETS: Add Forelimb Emergence Data (details.forelimb) to morph.data.mm -----------------------
# create column to store the developmental data
morph.data.mm$details.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.mm$unique.id.indiv))){
  morph.data.mm$details.forelimb[morph.data.mm$unique.id.indiv == unique(morph.data.mm$unique.id.indiv)[i]] <- devo.data$details.forelimb[devo.data$unique.id.indiv ==  unique(morph.data.mm$unique.id.indiv)[i]]
}

morph.data.mm$details.forelimb.binary = NA
morph.data.mm$details.forelimb.binary[morph.data.mm$details.forelimb == "left arm emerged"] = 0
morph.data.mm$details.forelimb.binary[morph.data.mm$details.forelimb == "right arm emerged"] = 1


# COMPILE DATASETS: Add Developmental Data (days to forelimb) and Mean Developmental Data to morph.data.mm -----------------------
# create column to store the developmental data
morph.data.mm$days.forelimb = NA
morph.data.mm$mean.days.forelimb = NA

# fill developmental data column - days.forelimb for each individual
for(i in 1:length(unique(morph.data.mm$unique.id.indiv))){
  morph.data.mm$days.forelimb[morph.data.mm$unique.id.indiv == unique(morph.data.mm$unique.id.indiv)[i]] <- devo.data$days.forelimb[devo.data$unique.id.indiv ==  unique(morph.data.mm$unique.id.indiv)[i]]
}

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.mm$unique.id))){
  morph.data.mm$mean.days.forelimb[morph.data.mm$unique.id == unique(morph.data.mm$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.mm$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.mm.juv -----------------------
# create column to store the developmental data
morph.data.mm.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.mm.juv$unique.id))){
  morph.data.mm.juv$mean.days.forelimb[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.mm.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.juv -----------------------
# create column to store the developmental data
morph.data.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to morph.data.mm -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.mm$mean.mm.mass.g = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.mm$unique.id[is.na(morph.data.mm$unique.id) == FALSE]))){
  morph.data.mm$mean.mm.mass.g[morph.data.mm$unique.id == unique(morph.data.mm$unique.id[is.na(morph.data.mm$unique.id) == FALSE])[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm$unique.id) == FALSE])[i]], na.rm = TRUE)
}

# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to morph.data.mm.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.mm.juv$mean.mm.mass.g = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE]))){
  morph.data.mm.juv$mean.mm.mass.g[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to morph.data.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.juv$mean.mm.mass.g = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.mm.mass.g[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Create scaled mass index body condition column  -----------------------
#he arithmetic mean of L is a suitable value for L0, and MiLi variables represent the raw data for each individual i.This index adjusts the mass of all individuals to that which they would have at length L0(represented by the vertical dashed line). 
# SMI = mass of individual(mean or median length of study population / length of individual) ^ scaling exponent
#first six for mm
plot(morph.data.mm.juv$mass.g~morph.data.mm.juv$svl.mm)

scalingeq = smatr::sma(log(morph.data.mm.juv$mass.g)~log(morph.data.mm.juv$svl.mm))
pearsons = cor.test(log(morph.data.mm.juv$mass.g), log(morph.data.mm.juv$svl.mm)) #unclear if pearsons should be calculated from log-transformed or untransformed data...it does change
scaling = scalingeq$coef[[1]][2,1]/pearsons$estimate
SUL.pop = mean(morph.data.mm.juv$svl.mm, na.rm = TRUE)

morph.data.mm$smi = morph.data.mm$mass.g*(SUL.pop/morph.data.mm$svl.mm)^scaling

morph.data.mm.juv$smi = morph.data.mm.juv$mass.g*(SUL.pop/morph.data.mm.juv$svl.mm)^scaling

morph.data.juv$smi = morph.data.juv$mass.g*(SUL.pop/morph.data.juv$svl.mm)^scaling


# COMPILE DATASETS: Add Average Metamorphosis body condition (mean smi at mm) to morph.data.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.juv$mean.mm.smi = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.mm.smi[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(morph.data.mm$smi[morph.data.mm$unique.id ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical smi data (i.e. poorer, average, better tanks) -- everyone seems to be better...maybe not useful at all
morph.data.juv$smi.cat = NA

for(i in 1:length(unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE]))){
  
  if(morph.data.juv$smi[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]][1] <= mean(morph.data.mm$smi, na.rm = TRUE) - (sd(morph.data.mm$smi, na.rm = TRUE)/2)){
    morph.data.juv$smi.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]] = "poorer"
  }
  
  if(morph.data.juv$smi[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]][1] >= mean(morph.data.mm$smi, na.rm = TRUE) + (sd(morph.data.mm$smi, na.rm = TRUE)/2)){
    morph.data.juv$smi.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]] = "better"
  }
}
morph.data.juv$smi.cat[is.na(morph.data.juv$smi.cat) == TRUE & is.na(morph.data.juv$juv.tank.id) == FALSE] = "average"
morph.data.juv$smi.cat = factor(morph.data.juv$smi.cat, levels = c("poorer", "average", "better"))


# COMPILE DATASETS: Add Average Metamorphosis body condition (mean smi at mm) to morph.data.mm -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.mm$mean.mm.smi = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.mm$mean.mm.smi[morph.data.mm$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(morph.data.mm$smi[morph.data.mm$unique.id ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Average Metamorphosis body condition (mean smi at mm) to morph.data.mm.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.mm.juv$mean.mm.smi = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.mm.juv$unique.id))){
  morph.data.mm.juv$mean.mm.smi[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]] <- mean(morph.data.mm$smi[morph.data.mm$unique.id ==  unique(morph.data.mm.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Average Gosner stage (avg.gosner.stage) to morph.data.tad -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.tad$avg.gosner.stage = NA

# fill avg.gosner.stage for morph.data.tad
for(i in 1:length(unique(morph.data.tad$unique.id))){
  morph.data.tad$avg.gosner.stage[morph.data.tad$unique.id == unique(morph.data.tad$unique.id)[i]] <- mean(devo.data.nonmm$avg.gosner.stage[devo.data.nonmm$unique.id ==  unique(morph.data.tad$unique.id)[i]], na.rm = TRUE)
}

# COMPILE DATASETS: Add Average Gosner stage (avg.gosner.stage) to morph.data.mm -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.mm$avg.gosner.stage = NA

# fill avg.gosner.stage for morph.data.tad
for(i in 1:length(unique(morph.data.mm$unique.id.larv))){
  morph.data.mm$avg.gosner.stage[morph.data.mm$unique.id.larv == unique(morph.data.mm$unique.id.larv)[i]] <- mean(devo.data.nonmm$avg.gosner.stage[devo.data.nonmm$unique.id ==  unique(morph.data.mm$unique.id.larv)[i]], na.rm = TRUE)
}


#COMPILE DATASETS: Create clutchtank for later analysis   -----------------------
morph.data.tad$clutchtank = paste(morph.data.tad$clutch, morph.data.tad$larv.tank.id, sep="_")
morph.data.mm$clutchtank = paste(morph.data.mm$clutch, morph.data.mm$larv.tank.id, sep="_")
morph.data.mm.juv$clutchtank = paste(morph.data.mm.juv$clutch, morph.data.mm.juv$juv.tank.id, sep="_")
morph.data.juv$clutchtank = paste(morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep="_")


# COMPILE DATASETS: Create tank-level means dataset from morph.data.juv and morph.data.mm.juv -----------------------
morph.data.juv.summ = morph.data.juv %>%
  group_by(post.mm.weeks, post.mm.weeks.num, gs, gs.code, clutch, treatment, juv.tank.id, clutchtank, mean.mm.mass.g, mean.days.forelimb, mean.mm.smi) %>%
  summarise(mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE),
            mean.r.forelimb.mm = mean(r.forelimb.mm, na.rm = TRUE),
            mean.r.tibia.mm = mean(r.tibia.mm, na.rm = TRUE),
            mean.r.thigh.mm = mean(r.thigh.mm, na.rm = TRUE),
            mean.smi = mean(smi, na.rm = TRUE)
  )

morph.data.mm.juv.summ = morph.data.mm.juv %>%
  group_by(post.mm.weeks, post.mm.weeks.num, gs, gs.code, clutch, treatment, juv.tank.id, clutchtank, mean.mm.mass.g, mean.days.forelimb, mean.mm.smi) %>%
  summarise(mean.mass.g = mean(mass.g, na.rm = TRUE),
            mean.svl.mm = mean(svl.mm, na.rm = TRUE),
            mean.r.forelimb.mm = mean(r.forelimb.mm, na.rm = TRUE),
            mean.r.tibia.mm = mean(r.tibia.mm, na.rm = TRUE),
            mean.r.thigh.mm = mean(r.thigh.mm, na.rm = TRUE),
            mean.smi = mean(smi, na.rm = TRUE)
  )



# COMPILE DATASETS: Create metrics vector to feed into plots and analyses later  -----------------------
metrics.tad = colnames(morph.data.tad)[17:22]
metrics.mm = colnames(morph.data.mm)[17:21]
metrics.mm.juv = colnames(morph.data.mm.juv)[c(11:15,21)]
metrics.juv = colnames(morph.data.juv)[c(24:28,33)]

#create summary table
morph.data.mm.summ = morph.data.mm %>% 
  group_by(treatment) %>% 
  summarise(min.days.forelimb = min(days.forelimb, na.rm = TRUE), 
            max.days.forelimb = max(days.forelimb, na.rm = TRUE),
            min.mass.g = min(mass.g, na.rm = TRUE), 
            max.mass.g = max(mass.g, na.rm = TRUE)
            )

morph.data.mm.summ$range.days.forelimb = morph.data.mm.summ$max.days.forelimb - morph.data.mm.summ$min.days.forelimb
morph.data.mm.summ$range.mass.g = morph.data.mm.summ$max.mass.g - morph.data.mm.summ$min.mass.g

#before analysis, remove from juvenile data anyone that didn't get placed in a juvenile tank. we want to still keep these in other dataframes, but don't want it here.
morph.data.juv = morph.data.juv[is.na(morph.data.juv$juv.tank.id) == FALSE,]

# PLOT DATASETS: Effect of rearing density on morphometrics during larval stage -----------------------

fig.3a.tl <- ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=tl.cm, x = week, color = treatment, fill =treatment), show.legend = F) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 1, color = "black", aes(x = week, y = tl.cm, group = treatment), show.legend = F) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = week, y = tl.cm, group = treatment), show.legend = F) +
  
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = tl.cm, fill=treatment), show.legend = T) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "total length (cm)", limits = c(1.4,4.5), breaks = seq(1.5,4.5,0.5)) +
  scale_x_continuous(name = "larval age (weeks)", labels = seq(0,10,by=1), breaks = seq(0, 10, by = 1))


fig.3b.bl <- ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=bl.cm, x = week, color = treatment, fill =treatment), show.legend = F) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 1, color = "black", aes(x = week, y = bl.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = week, y = bl.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = bl.cm, fill=treatment), show.legend = F) + 

  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body length (cm)", limits = c(0.3,1.7), breaks = seq(0.4,1.7,0.3)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


fig.3c.tail <- ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=tail.cm, x = week, color = treatment, fill =treatment), show.legend = F) +

  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = tail.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = week, y = tail.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = tail.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "tail length (cm)", limits = c(0.4,3.3), breaks = seq(0.6,3.3,0.6)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


# temporary plot to assess proportion of total length that comprises body length
ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=bl.cm/tl.cm, x = week, color = treatment, fill =treatment), show.legend = F) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = bl.cm/tl.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = week, y = bl.cm/tl.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = bl.cm/tl.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=19.5, color = "black"), 
        axis.text.y=element_text(size=19.5, color = "black"), 
        axis.title.x = element_text(size=19.5),
        axis.title.y = element_text(size=19.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body length/total length") + 
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


# temporary plot to assess proportion of total length that comprises tail length
ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=tail.cm/tl.cm, x = week, color = treatment, fill =treatment), show.legend = F) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = tail.cm/tl.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 0.8, colour="black", alpha=1, aes(x = week, y = tail.cm/tl.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = tail.cm/tl.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "tail length/total length") + 
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))




fig.3a.hw <- ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=hw.cm, x = week, color = treatment, fill =treatment), show.legend = F) +
  
  #tank means
  stat_summary(data = morph.data.tad, fun = mean, geom="point", pch=21, size=3, aes(x = week, y = hw.cm, fill=treatment, group = clutchtank, color = treatment), show.legend = F) + 
  stat_summary(data = morph.data.tad, fun = mean, geom="line", pch=21, size=0.5, aes(x = week, y = hw.cm, group = clutchtank, color = treatment), show.legend = F) + 
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = hw.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(x = week, y = hw.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = hw.cm, fill=treatment), show.legend = F) + 

  
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "head width (cm)")  + #, limits = c(1,4.5), breaks = seq(1,4.5,0.5)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=bw.cm, x = week, color = treatment, fill =treatment), show.legend = T) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = bw.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(x = week, y = bw.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = bw.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body width")  + #, limits = c(1,4.5), breaks = seq(1,4.5,0.5)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


fig.3a.bw.bl <- ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=bw.cm/bl.cm, x = week, color = treatment, fill =treatment), show.legend = T) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = bw.cm/bl.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(x = week, y = bw.cm/bl.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = bw.cm/bl.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "body width/body length")  + #, limits = c(1,4.5), breaks = seq(1,4.5,0.5)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))


ggplot() + 
  
  #individual tanks
  geom_jitter(size = 1, alpha = 0.7, data = morph.data.tad, pch = 21, aes(y=tmw.cm, x = week, color = treatment, fill =treatment), show.legend = T) +
  
  #treatment means
  stat_summary(data = morph.data.tad, fun=mean, geom="line", size = 0.8, color = "black", aes(x = week, y = tmw.cm, group = treatment), show.legend = FALSE) + #show.legend is FALSE so that alpha = 1 legend from stat summary will be plotted
  stat_summary(data = morph.data.tad, fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(x = week, y = tmw.cm, group = treatment), show.legend = F) +
  stat_summary(data = morph.data.tad, fun = mean, geom="point", color = "black", pch=21, size=5, stroke = 1, aes(x = week, y = tmw.cm, fill=treatment), show.legend = F) + 
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "tail-muscle width")  + #, limits = c(1,4.5), breaks = seq(1,4.5,0.5)) +
  scale_x_continuous(name = "larval age (weeks)", breaks = seq(0, 12, by = 1))




# PLOT DATASETS: Effect of rearing density on morphometrics AT metamorphosis for first six individuals from each tank -----------------------
fig.4a <- ggplot(data = morph.data.mm, aes(y=mass.g, x = treatment, fill = treatment, color = treatment)) + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), size = 1, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(alpha = 0.6, size = 1, show.legend = FALSE, aes(fill = treatment), outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  geom_signif(comparisons = list(c("low density", "high density")), 
              map_signif_level=TRUE,
              annotations = "ns",
              size = 1,
              tip_length = 0,
              vjust = -0.3,
              textsize = 6,
              color = "black") + 
  
  scale_y_continuous(name = "metamorphic mass (g)", limits = c(0.07,0.29), breaks = seq(0.08,0.28,0.05)) +
  scale_x_discrete(name = "")

#fig.2b <- ggExtra::ggMarginal(fig.2b, margins = "y", groupColour = TRUE, groupFill = TRUE)

fig.4b <- ggplot(data = morph.data.mm, aes(y=smi, x = treatment, fill = treatment, color = treatment)) + 
 
   geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), size = 1, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(alpha = 0.6, size = 1, show.legend = FALSE, aes(fill = treatment), outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  geom_signif(comparisons = list(c("low density", "high density")), 
              map_signif_level=TRUE,
              annotations = "ns",
              size = 1,
              tip_length = 0,
              vjust = -0.3,
              textsize = 6,
              y_position = 0.26, 
              color = "black") + 
  
  scale_y_continuous(name = "metamorphic SMI", limits = c(0.135,0.2695), breaks = seq(0.11,0.28,0.03)) +
  scale_x_discrete(name = "")

#fig.2c <- ggExtra::ggMarginal(fig.2c, margins = "y", groupColour = TRUE, groupFill = TRUE)


# PLOT DATASETS: Effect of rearing density on morphometrics at and after metamorphosis for first six individuals from each tank -----------------------
# option 1 = x-y plot with summarized mean and +/- 1 se for all metrics

# create vector of morphometrics and y-axis labels
yaxis.names = c("mass (g)", "snout-vent length (mm)", "forearm length (mm)", "tibia length (mm)", "thigh length (mm)", "scaled mass index")

# create empty list to fill with morphometrics plot
plotList.morph = vector(mode = "list", length = length(metrics.mm.juv))

# fill plot list with each metric for mass and svl
for(i in 1:length(metrics.mm.juv[c(1,2)])){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics.mm.juv[i]]], x = post.mm.weeks, color = treatment)) + 
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
for(i in 3:length(metrics.mm.juv)){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=.data[[metrics.mm.juv[i]]], x = post.mm.weeks, color = treatment)) + 
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
    scale_y_continuous(name = yaxis.names[i]) +
    scale_x_discrete(name = "post-metamorphic age (weeks)")
}

# add the legend as the final plot within plot list so that it can be graphed within the grid
#plotList.morph[[length(plotList.morph) + 1]] <- as_ggplot(get_legend(plotList.morph))

# create panel plot with only mass and svl data across all sampling points
ggarrange(plotlist = plotList.morph, 
          common.legend = FALSE,
          legend = "bottom",
          labels = c("a", "b", "c", "d", "e", "f"),
          font.label = list(size = 20, color = "black"))


plot.morph2 <- ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=mass.g, x = post.mm.weeks, color = treatment)) + 
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
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mass (g)") +
  scale_x_discrete(name = "post-metamorphic age (weeks)")


# create panel plot with all other morphometrics data (not mass or svl) across all sampling points
ggarrange(plotlist = plotList.morph[c(3,4,5)], 
          ncol = 1,
          nrow = 3,
          common.legend = TRUE,
          legend = "right",
          labels = c("a", "b", "c", "d", "e", ""),
          font.label = list(size = 16, color = "black"))


# PLOT DATASETS: Effect of rearing density on correlation between morphology metrics across larval duration -----------------------

# create empty list to fill with morphometrics plot
plotList.morph = vector(mode = "list", length = length(metrics.tad))
yaxis.names.tad = c("total length (cm)", "body length (cm)", "tail length (cm)", "head width (cm)", "body width (cm)", "tail-muscle width (cm)")

# fill plot list with each metric
for(i in 1:length(metrics.tad)){
  plotList.morph[[i]] <- 
    ggplot(data = morph.data.tad[!is.na(morph.data.tad$larv.tank.id),], aes(y=.data[[metrics.tad[i]]], x = tl.cm, color = treatment)) + 
    facet_grid(rows = vars(factor(week))) +
    geom_abline(slope = 1, intercept = 0, na.rm = FALSE, show.legend = NA, linetype = 2) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
    scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
    scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=20),
          axis.text.x=element_text(size=12, color = "black"), 
          axis.text.y=element_text(size=12, color = "black"), 
          axis.title.x=element_text(size=12, color = "black"), 
          axis.title.y = element_text(size=12),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    #scale_y_continuous(name = yaxis.names.tad[i], limits = c(0,4)) +
    scale_x_continuous(name = "total length (cm)")
}

# create panel plot with all morphometrics data across all sampling points
ggarrange(plotlist = plotList.morph,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b", "c", "d", "e"),
          font.label = list(size = 20, color = "black"))

plotList.morph


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


# PLOT DATASETS: Correlation matrix between morphology metrics during larval period -----------------------

ggplot(data = morph.data.corr.tad, aes(y=estimate, x = factor(week), color = metric)) + 
  facet_grid(rows = vars(treatment)) +
  geom_point(size = 3, alpha = 1) +
  geom_line(size = 1, alpha = 1, aes(y=estimate, x = factor(week), color = metric, group = metric)) +
  # scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  # scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
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
  scale_x_discrete(name = "metric")


# PLOT DATASETS: Correlation matrix between morphology metrics at and after metamorphosis -----------------------

morph.data.corr$metric = factor(morph.data.corr$metric, levels = c("svl.mm", "r.forelimb.mm", "r.tibia.mm", "r.thigh.mm"))

ggplot(data = morph.data.corr, aes(y=estimate, x = post.mm.weeks, color = metric)) + 
  #facet_grid(rows = vars(treatment)) +
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
  scale_x_discrete(name = "post-metamorphic age")



# PLOT DATASETS: Effect of developmental speed and treatment on morphometrics & SMI AT metamorphosis ---------------------

fig.5a <- ggplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, seed = 0), pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm,
             colour = "black", stroke = 0.75,
             aes(y=mass.g, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  geom_ribbon(data = predicted.df.morph.mm,
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.morph.mm, size = 1, color = "black", aes(x = x, y = predicted)) +
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = -0.001", sep = ""))),
           x=43.5, y = 0.28, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 11.4", sep = ""))),
           x=43.5, y = 0.269, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.001",
           x=43.5, y = 0.254, color = "black", size = 5, hjust=0) +

  scale_fill_manual(values=c("gray65", "gray65")) +
  
  theme_bw() +
  theme(legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic mass (g)", limits = c(0.08,0.28), breaks = seq(0.08,0.28,0.05)) +
  scale_x_continuous(name = "larval duration (days)", limits = c(43.5, 73.5), breaks = seq(44, 75, by = 5))


fig.5b <- ggplot() + 
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, seed = 0), pch = 21, alpha = 0.7, size = 3, show.legend = FALSE,
             data = morph.data.mm,
             color = "black", stroke = 0.75,
             aes(y=smi, x = days.forelimb, fill = treatment)) +
  geom_ribbon(data = predicted.df.smi.mm,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.smi.mm, size = 1, color = "black", aes(x = x, y = predicted)) +
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = +0.0005", sep = ""))),
           x=43.5, y = 0.26, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 6.0", sep = ""))),
           x=43.5, y = 0.253, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.05",
           x=43.5, y = 0.246, color = "black", size = 5, hjust=0) +
 
  scale_fill_manual(values=c("gray65", "gray65")) +

  theme_bw() +
  theme(legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic SMI", limits = c(0.14,0.26), breaks = seq(0.11,0.28,0.03)) +
  scale_x_continuous(name = "larval duration (days)", limits = c(43.5, 73.5), breaks = seq(44, 75, by = 5))



fig.5c <- ggplot() + 
  geom_point(pch = 21, alpha = 0.7, show.legend = FALSE, size = 3,
             color = "black", stroke = 0.75, 
             data = morph.data.mm,
             aes(y=smi, x = mass.g, fill = treatment)) +
  geom_ribbon(data = predicted.df.smi.mm.2,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.smi.mm.2, size = 1, color = "black", aes(x = x, y = predicted)) +
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = +0.119", sep = ""))),
           x=0.08, y = 0.26, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 4.2", sep = ""))),
           x=0.08, y = 0.253, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.001",
           x=0.08, y = 0.246, color = "black", size = 5, hjust=0) +
  
  scale_fill_manual(values=c("gray65", "gray65")) +
  
  theme_bw() +
  theme(axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic SMI", limits = c(0.14,0.26), breaks = seq(0.14, 0.26, 0.03)) +
  scale_x_continuous(name = "metamorphic mass (g)", limits = c(0.08,0.28), breaks = seq(0.08,0.28,0.05))



# PLOT DATASETS: Effect of metamorphic mass on morphometrics at and after metamorphosis ---------------------

#plotting individual mass on y-axis but no means; most similar to model - response as individual but no means shown
fig.6a <- ggplot() +
  
  geom_point(data = morph.data.mm.juv[is.na(morph.data.mm.juv$juv.tank.id) == FALSE,], #do not plot values for which juvenile tank doesn't exist
             position = position_jitter(width = 0.7, height = NULL, seed = 0.1),
             pch = 21, size = 2.5, color = "black", alpha = 0.75,
             aes(x = post.mm.weeks.num, y = mass.g, group = clutchtank, fill = mean.mm.mass.g), show.legend = T) +
  
  scale_color_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  scale_fill_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  
  labs(fill = "tank mean metamorphic mass (g)") +
  
  #add statistical results
  annotate(geom = "text", 
           label = as.character(expression(paste(beta, " = +4.97", sep = ""))),
           x=30, y = 0.3, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = as.character(expression(paste(chi^2, " = 6.62", sep = ""))),
           x=30, y = 0.2, color = "black", size = 5, hjust=0, parse = TRUE) +
  annotate(geom = "text", 
           label = "p < 0.05",
           x=30, y = 0.1, color = "black", size = 5, hjust=0) +

theme_bw() +
  theme(legend.position = c(0.02,0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "individual mass (g)", limits = c(0.05,1.1), breaks = seq(0.05, 1.1, 0.25)) + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))

#plotting tank-level means on y-axis
fig.6b <- ggplot() +

  stat_summary(data = morph.data.mm.juv[is.na(morph.data.mm.juv$juv.tank.id) == FALSE,], fun = mean, #do not plot values for which juvenile tank doesn't exist 
               geom = "errorbar", position = position_jitter(width = 0.7, height = NULL, seed = 0.1),
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.5, alpha=0.7, aes(x = post.mm.weeks.num, y = mass.g, group = clutchtank), show.legend = F) +

  stat_summary(data = morph.data.mm.juv[is.na(morph.data.mm.juv$juv.tank.id) == FALSE,], fun = mean,
               geom = "line", size = 0.75, position = position_jitter(width = 0.7, height = NULL, seed = 0.1),
               alpha=0.75, aes(x = post.mm.weeks.num, y = mass.g, group = clutchtank, color = mean.mm.mass.g), show.legend = F) +
  
  stat_summary(data = morph.data.mm.juv[is.na(morph.data.mm.juv$juv.tank.id) == FALSE,], fun = mean, 
               geom = "point", pch = 21, position = position_jitter(width = 0.7, height = NULL, seed = 0.1), size = 2.5, color = "black", stroke = 0.5,
               alpha=0.75, aes(x = post.mm.weeks.num, y = mass.g, group = clutchtank, fill = mean.mm.mass.g), show.legend = T) +

  scale_color_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  scale_fill_gradient2(low = alpha("#FFF7EC", 0.9), mid = alpha("#FC8D59", 1), high = alpha("darkred", 1), midpoint = 0.161875) +
  
  labs(fill = "tank mean metamorphic mass (g)") +
  
  theme_bw() +
  theme(legend.position = c(0.02,0.99),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "tank mean mass (g)", limits = c(0.05,1.1), breaks = seq(0.05, 1.1, 0.25)) + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))



fig.XXX <- ggplot() + 
  
  geom_point(data = morph.data.mm.juv[morph.data.mm.juv$post.mm.weeks.num == 1.5,], 
             pch = 21, size = 3, color = "black", stroke = 0.5,
               alpha=1, aes(x = mean.mm.mass.g, y = mass.g, group = clutchtank, fill = mean.mm.mass.g), show.legend = F) + 
  
  scale_fill_gradient(low = "gray90", high = "gray1") +
  labs(fill = "tank-level mean metamorphic mass (g)") +

  #add descriptive text
  annotate(geom = "text", 
           label = "1-2 weeks",
           x=0.13, y = 0.25, color = "black", size = 5, hjust=0) +

  theme_bw() +
  theme(legend.position = c(0.02,0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "mass (g)") +
  scale_x_continuous(name = "tank-level mean metamorphic mass (g)", limits = c(0.13,0.195), breaks = seq(0.13,0.2,0.02), labels = seq(0.13,0.2,0.02))

fig.5c <- ggplot() + 
  
  geom_point(data = morph.data.mm.juv[morph.data.mm.juv$post.mm.weeks.num == 5,], 
             pch = 21, size = 3, color = "black", stroke = 0.5,
             alpha=1, aes(x = mean.mm.mass.g, y = mass.g, group = clutchtank, fill = mean.mm.mass.g), show.legend = F) + 
  
  scale_fill_gradient(low = "gray90", high = "gray1") +
  labs(fill = "tank-level mean metamorphic mass (g)") +
  
  #add descriptive text
  annotate(geom = "text", 
           label = "4-6 weeks",
           x=0.13, y = 0.35, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.position = c(0.02,0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "mass (g)") +
  scale_x_continuous(name = "tank-level mean metamorphic mass (g)", limits = c(0.13,0.195), breaks = seq(0.13,0.2,0.02), labels = seq(0.13,0.2,0.02))

fig.5d <- ggplot() + 
  
  geom_point(data = morph.data.mm.juv[morph.data.mm.juv$post.mm.weeks.num == 30.5,], 
             pch = 21, size = 3, color = "black", stroke = 0.5,
             alpha=1, aes(x = mean.mm.mass.g, y = mass.g, group = clutchtank, fill = mean.mm.mass.g), show.legend = F) + 
  
  scale_fill_gradient(low = "gray90", high = "gray1") +
  labs(fill = "tank-level mean metamorphic mass (g)") +
  
  #add descriptive text
  annotate(geom = "text", 
           label = "30-31 weeks",
           x=0.13, y = 1.65, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.position = c(0.02,0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "mass (g)") + 
  scale_x_continuous(name = "tank-level mean metamorphic mass (g)", limits = c(0.13,0.195), breaks = seq(0.13,0.2,0.02), labels = seq(0.13,0.2,0.02))







# PLOT DATASETS: Effect of metamorphic body condition on body condition at and after metamorphosis ---------------------
ggplot(data=morph.data.mm.juv,
       aes(y=smi, x = mass.g, color = mean.mm.smi, group = clutchtank)) +
  geom_point() +
  facet_grid(rows = vars(post.mm.weeks)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "lm") +
  geom_smooth(se = T, size = 0.8, alpha = 0.6, method = "lm", color = "red", data=morph.data.mm.juv,
              aes(y=smi, x = mass.g, color = mean.mm.smi, group = treatment, linetype = treatment))

# assess tank-level vs. individual-level fits with mean smi
plot.smi.1 <- ggplot() + 
  geom_jitter(width = 0.8, size = 2, alpha = 1, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.smi, group = clutchtank)) + #individual level response colored by mean tank mass at metamorphosis
  # scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  geom_line(data = morph.data.mm.juv.summ, pch = 21, size = 0.8, aes(y=mean.smi, x = post.mm.weeks.num, color = mean.mm.smi, group = clutchtank)) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_fill_gradient(low = "gray85", high = "gray1") +
  
  new_scale_fill() +
  # geom_ribbon(data = predicted.df.bodycond.juv,
  #             alpha = 0.5,
  #             mapping = aes(x = x, y = predicted, fill = group, group = group, ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_point(data = morph.data.mm.juv.summ, pch = 21, alpha = 1, size = 4, stroke = 1, aes(y=mean.smi, x = post.mm.weeks.num, fill = treatment, color = mean.mm.smi,  group = clutchtank)) + #tank level response
  
  new_scale_color() +
  geom_line(data = predicted.df.bodycond.juv, linetype = 2, size = 2, aes(x = x, y = predicted, color = group)) + #from model on individual smi
  geom_line(data = predicted.df.bodycond.juv.summ, color = "black", size = 2, linetype = 3, aes(x = x, y = predicted)) + #from model on tank-level means for smi
  
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  theme(legend.justification = c("left", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))

# assess tank-level vs. individual-level fits with mean days forelimb
plot.smi.2 <-ggplot() + 
  geom_jitter(width = 0.8, size = 2, alpha = 1, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.days.forelimb, group = clutchtank)) + #individual level response colored by mean tank mass at metamorphosis
  # scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  geom_line(data = morph.data.mm.juv.summ, pch = 21, size = 0.8, aes(y=mean.smi, x = post.mm.weeks.num, color = mean.days.forelimb, group = clutchtank)) +
  scale_color_gradient(low = "gray85", high = "gray1", limits = c(47,71)) +
  scale_fill_gradient(low = "gray85", high = "gray1", limits = c(47,71)) +
  
  new_scale_fill() +
  # geom_ribbon(data = predicted.df.bodycond.juv,
  #             alpha = 0.5,
  #             mapping = aes(x = x, y = predicted, fill = group, group = group, ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_point(data = morph.data.mm.juv.summ, pch = 21, alpha = 1, size = 4, stroke = 1, aes(y=mean.smi, x = post.mm.weeks.num, fill = treatment, color = mean.days.forelimb,  group = clutchtank)) + #tank level response
  
  new_scale_color() +
  geom_line(data = predicted.df.bodycond.juv, linetype = 2, size = 2, aes(x = x, y = predicted, color = group)) + #from model on individual smi
  geom_line(data = predicted.df.bodycond.juv.summ, color = "black", size = 2, linetype = 3, aes(x = x, y = predicted)) + #from model on tank-level means for smi
  
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  theme(legend.justification = c("left", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))



# x-y plot with observed and glmm-predicted values for mass
plot.smi.1 <- ggplot() + 
  facet_grid(cols=vars(treatment)) +
  # geom_ribbon(data = predicted.df.juv,
  #             alpha = 0.75,
  #             mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("large", "med", "small")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.smi, group = clutchtank)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "lm", data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.smi, group=clutchtank)) +
  #geom_line(data = predicted.df.juv, size = 1, aes(x = x.num, y = predicted, linetype=factor(group, levels = c("large", "med", "small")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(color="mean tank smi") +
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
  scale_y_continuous(name = "individual smi") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))

plot.smi.2 <- ggplot() + 
  facet_grid(cols=vars(treatment)) +
  # geom_ribbon(data = predicted.df.juv,
  #             alpha = 0.75,
  #             mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("large", "med", "small")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.mass.g, group = clutchtank)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "lm", data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.mass.g, group=clutchtank)) +
  #geom_line(data = predicted.df.juv, size = 1, aes(x = x.num, y = predicted, linetype=factor(group, levels = c("large", "med", "small")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(color="mean metamorphic mass (g)") +
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
  scale_y_continuous(name = "individual smi") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))

plot.smi.3 <- ggplot() + 
  facet_grid(cols=vars(treatment)) +
  # geom_ribbon(data = predicted.df.juv,
  #             alpha = 0.75,
  #             mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("large", "med", "small")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.days.forelimb, group = clutchtank)) +
  geom_smooth(se = F, size = 0.8, alpha = 0.7, method = "lm", data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.days.forelimb, group=clutchtank)) +
  #geom_line(data = predicted.df.juv, size = 1, aes(x = x.num, y = predicted, linetype=factor(group, levels = c("large", "med", "small")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(color="mean larval duration (days)") +
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
  scale_y_continuous(name = "individual smi") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))

ggarrange(plot.smi.1, plot.smi.2, plot.smi.3,
          ncol = 1,
          nrow = 3,
          common.legend = FALSE)


ggplot() + 
  geom_ribbon(data = predicted.df.bodycond.juv,
              alpha = 0.75,
              mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("better", "average", "poorer")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks.num, color = mean.mm.smi)) +
  geom_line(data = predicted.df.bodycond.juv, size = 1, aes(x = as.numeric(x.num), y = predicted, linetype=factor(group, levels = c("better", "average", "poorer")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(linetype="mean tank smi category", color="mean tank smi") +
  theme(#legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))


ggplot() + 
  # geom_ribbon(data = predicted.df.bodycond.juv,
  #             alpha = 0.75,
  #             mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("better", "average", "poorer")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  facet_grid(rows = vars(treatment)) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=smi, x = post.mm.weeks, color = mean.mm.smi)) +
  geom_line(data = predicted.df.bodycond.juv, size = 1, aes(x = as.numeric(x.num), y = predicted, linetype=factor(group, levels = c("better", "average", "poorer")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(linetype="mean tank smi category", color="mean tank smi") +
  theme(#legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.title = element_text(size=10),
    legend.text = element_text(size=10),
    axis.text.x=element_text(size=12, color = "black"), 
    axis.text.y=element_text(size=12, color = "black"), 
    axis.title.x=element_text(size=12, color = "black"), 
    axis.title.y = element_text(size=12),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))


# growth rate
ggplot() + 
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = growth.data.mm.juv, aes(y=growth.mass.g, x = post.mm.weeks.num, color = mean.mm.mass.g)) +
  geom_ribbon(data = predicted.df.growth.juv,
              alpha = 0.75,
              mapping = aes(x = x.num, y = predicted, ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.growth.juv, size = 1, aes(x = x.num, y = predicted)) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  theme_bw() +
  labs(color="mean tank metamorphic mass (g)") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "growth rate (g/day)") +
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))



# PLOT DATASETS: Effect of larval duration on morphometrics at and after metamorphosis ---------------------

# x-y plot with observed and glmm-predicted values for mass
plot.morph4 <- ggplot() + 
  #geom_ribbon(data = predicted.df.juv,
              # alpha = 0.75,
              # mapping = aes(x = as.numeric(x.num), y = predicted, fill=factor(group, levels = c("large", "med", "small")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2.5, alpha = 0.8, data = morph.data.mm.juv, aes(y=mass.g, x = post.mm.weeks.num, color = mean.days.forelimb)) +
  #geom_line(data = predicted.df.juv, size = 1, aes(x = x.num, y = predicted, linetype=factor(group, levels = c("large", "med", "small")))) +
  scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
  labs(linetype="mean tank metamorphic mass (category)", color="mean larval duration (days)") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "individual metamorphic mass (g)") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))


# PANEL PLOTS + EXPORT PLOTS: Create panel plot with mass as a function of treatment and larval duration across all sampling points --------------
# pdf("~/Desktop/R Working Directory/Plots/Figure3.pdf", width = 11, height = 8.5)
# ggarrange(fig.3a.tl, 
#           ggarrange(fig.3b.bl, fig.3c.tail,
#                     ncol = 1,
#                     nrow = 2,
#                     labels = c("", "c"),
#                     font.label = list(size = 20, color = "black")),
#           ncol = 2,
#           nrow = 1,
#           common.legend = TRUE,
#           legend = "bottom",
#           labels = c("a", "b"),
#           font.label = list(size = 20, color = "black"))
# dev.off()

png("~/Desktop/R Working Directory/Plots/Figure3.png", units = "in", res = 300, width = 8.5, height = 8.5)
ggarrange(fig.3a.tl, 
          ggarrange(fig.3b.bl, fig.3c.tail,
                    ncol = 2,
                    nrow = 1,
                    labels = c("", "c"),
                    font.label = list(size = 20, color = "black")),
          ncol = 1,
          nrow = 2,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()

png("~/Desktop/R Working Directory/Plots/Figure4.png", units = "in", res = 300, width = 8.5, height = 11)
ggarrange(fig.4a, fig.4b,
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()

png("~/Desktop/R Working Directory/Plots/Figure5.png", units = "in", res = 300, width = 11, height = 8.5)
ggarrange(fig.5a,
          ggarrange(fig.5b, fig.5c,
                    ncol = 2,
                    nrow = 1,
                    labels = c("b", "c"),
                    font.label = list(size = 20, color = "black")),
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          labels = c("a", ""),
          font.label = list(size = 20, color = "black"))
dev.off()

png("~/Desktop/R Working Directory/Plots/Figure6.png", units = "in", res = 300, width = 16, height = 8.5)
ggarrange(fig.6a, fig.6b,
          ncol = 1,
          nrow = 2,
          common.legend = FALSE,
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))
dev.off()



# DOWNLOAD DATA FOR JOURNAL SUBMISSION  ---------------------
#set wd
setwd("~/Desktop/R Working Directory/2023_ranid_devo_plasticity/Submission_RoySocB")

write.csv(morph.data.tad[,c(1:2, 5:11, 17:22, 33, 35)], file = "Database_Morphometrics_Tadpole.csv", row.names = FALSE)

#morph.data.juv$post.mm.weeks = as.character(morph.data.juv$post.mm.weeks)
write.csv(morph.data.juv[,c(1,3,10:15,17, 24:28,31,32:36, 38)], file = "Database_Morphometrics_Juvenile.csv", row.names = FALSE)

write.csv(morph.data.mm[,c(1,6:13, 17:21, 27:29, 32:36, 38)], file = "Database_Morphometrics_Metamorph.csv", row.names = FALSE)

write.csv(morph.data.mm.juv[,c(1, 3, 5:15, 17, 19:23)], file = "Database_Morphometrics_MetamorphAndJuvenile.csv", row.names = FALSE)


# SUMMARY TABLES: Morphological metrics at and after metamorphosis ---------------------


summ.morph.data.tad = morph.data.tad %>%
  group_by(week, treatment) %>%
  summarise(mean.tmw.cm = round(mean(tmw.cm, na.rm = TRUE), 4),
            sd.tmw.cm  = round(sd(tmw.cm, na.rm = TRUE), 3),
            n = n()
  )

#reset working directory to be outputs folder
setwd("~/Desktop/R Working Directory/Outputs")

# create summary table for morphometrics that includes rearing density for first six
morph.data.mm.summ = morph.data.mm %>%
  group_by(treatment) %>%
  summarise(mean.mass = round(mean(mass.g, na.rm = TRUE), 4),
            sd.mass = round(sd(mass.g, na.rm = TRUE), 3),
            mean.smi = round(mean(smi, na.rm = TRUE), 4),
            sd.smi= round(sd(smi, na.rm = TRUE), 3),
            mean.svl = round(mean(svl.mm, na.rm = TRUE), 2),
            sd.svl = round(sd(svl.mm, na.rm = TRUE), 2),
            mean.forearm = round(mean(r.forelimb.mm, na.rm = TRUE), 2),
            sd.forearm = round(sd(r.forelimb.mm, na.rm = TRUE), 2),
            mean.tibia = round(mean(r.tibia.mm, na.rm = TRUE), 2),
            sd.tibia = round(sd(r.tibia.mm, na.rm = TRUE),2),
            mean.thigh = round(mean(r.thigh.mm, na.rm = TRUE),2),
            sd.thigh = round(sd(r.thigh.mm, na.rm = TRUE),2)
  )

# create summary table for morphometrics that includes water reduction for first six
morph.data.mm.water.reduc.summ = morph.data.mm %>%
  group_by(water.level.reduc) %>%
  summarise(mean.mass = round(mean(mass.g, na.rm = TRUE), 4),
            sd.mass = round(sd(mass.g, na.rm = TRUE), 3),
            mean.smi = round(mean(smi, na.rm = TRUE), 4),
            sd.smi= round(sd(smi, na.rm = TRUE), 3),
            mean.svl = round(mean(svl.mm, na.rm = TRUE), 2),
            sd.svl = round(sd(svl.mm, na.rm = TRUE), 2),
            mean.forearm = round(mean(r.forelimb.mm, na.rm = TRUE), 2),
            sd.forearm = round(sd(r.forelimb.mm, na.rm = TRUE), 2),
            mean.tibia = round(mean(r.tibia.mm, na.rm = TRUE), 2),
            sd.tibia = round(sd(r.tibia.mm, na.rm = TRUE),2),
            mean.thigh = round(mean(r.thigh.mm, na.rm = TRUE),2),
            sd.thigh = round(sd(r.thigh.mm, na.rm = TRUE),2)
  )


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


# create summary table for forelimb emergence that doesn't include rearing density for first six
morph.data.mm.forelimb = morph.data.mm %>%
  group_by(details.forelimb) %>%
  summarise(n = n(),
)

# create summary table for forelimb emergence that doesn't include rearing density for all
morph.data.mm.all = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica and remove NAs
morph.data.mm.all = morph.data.mm.all[morph.data.mm.all$gs.code == "RS" & is.na(morph.data.mm.all$gs.code) == FALSE,]

# subset devo dataset to only include first six individuals and non-overflow individuals
morph.data.mm.all = morph.data.mm.all[morph.data.mm.all$treatment != "overflow",]

morph.data.mm.forelimb.all = morph.data.mm.all %>%
  group_by(details.forelimb) %>%
  summarise(n = n(),
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


# create summary table for sampling
sampling.summ = morph.data.mm %>%
  group_by(treatment, clutch) %>%
  summarise(replicates = length(unique(larv.tank.id)))
sampling.summ$ind.per.rep[sampling.summ$treatment == "low density"] = 40
sampling.summ$ind.per.rep[sampling.summ$treatment == "high density"] = 100
sampling.summ$ind.tot = sampling.summ$replicates*sampling.summ$ind.per.rep











# TEMPORARY: Figure out how many completed measurements we have for each tank across weeks  -----------------------
morph.data.summ = morph.data.tad %>%
  group_by(gs.code, clutch, treatment, larv.tank.id) %>%
  filter(is.na(tl.cm) == FALSE) %>%
  summarise(unique.weeks = length(unique(week)),
            weeks = paste(unique(week), collapse=","))


# TEMPORARY: Figure out relationship between mean larval duration, mean metamorphic mass, and mean metamorphic smi  -----------------------
# as larval duration increases, mass decreases, body condition increases
ggarrange(ggplot(data=morph.data.mm, aes(x=mean.days.forelimb, y=mean.mm.mass.g, color = treatment)) +
  geom_point() +
  #facet_grid(cols=vars(treatment)) +
  stat_smooth(method = "lm"),
  
  ggplot(data=morph.data.mm, aes(x=days.forelimb, y=mass.g, color = treatment)) +
    geom_point() +
    #facet_grid(cols=vars(treatment)) +
    stat_smooth(method = "lm"),
  
  ncol = 2, nrow = 1)

ggarrange(ggplot(data=morph.data.mm, aes(x=mean.mm.mass.g, y=mean.mm.smi, color = treatment)) +
  geom_point() +
  #facet_grid(cols=vars(treatment)) +
  stat_smooth(method = "lm"),
  
  ggplot(data=morph.data.mm, aes(x=mass.g, y=smi, color = treatment)) +
    geom_point() +
    #facet_grid(cols=vars(treatment)) +
    stat_smooth(method = "lm"),
  
  ncol = 2, nrow = 1)


ggarrange(ggplot(data=morph.data.mm, aes(x=mean.days.forelimb, y=mean.mm.smi, color = treatment)) +
  geom_point() +
  #facet_grid(cols=vars(treatment)) +
  stat_smooth(method = "lm"),
  
  ggplot(data=morph.data.mm, aes(x=days.forelimb, y=smi, color = treatment)) +
    geom_point() +
    #facet_grid(cols=vars(treatment)) +
    stat_smooth(method = "lm"),
  
  ncol = 2, nrow = 1)


ggarrange(ggplot(data=morph.data.mm, aes(x=clutchtank, y=mass.g, color = treatment)) +
            facet_grid(cols = vars(treatment)) +
            geom_boxplot() +
            scale_y_continuous(limits = c(0.05,0.3)),
          
          ggplot(data=morph.data.mm, aes(x=clutchtank, y=mean.mm.mass.g, color = treatment)) +
            facet_grid(cols = vars(treatment)) +
            geom_boxplot() +
            scale_y_continuous(limits = c(0.05,0.3)),
          
          ncol = 2, nrow = 1)

ggarrange(ggplot(data=morph.data.mm, aes(x=days.forelimb, y=smi, color = treatment)) +
  facet_grid(cols = vars(treatment)) +
  geom_point() +
    stat_smooth(method = "lm") +
    scale_y_continuous(limits = c(0.10, 0.25)),
  
ggplot(data=morph.data.mm, aes(x=mean.days.forelimb, y=mean.mm.smi, color = treatment)) +
  facet_grid(cols = vars(treatment)) +
  geom_point() +
  stat_smooth(method = "lm") +
  scale_y_continuous(limits = c(0.10, 0.25)),

ncol = 2, nrow = 1)
  

summary(lmer(mean.mm.smi ~ mean.days.forelimb + (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], REML = FALSE))
summary(lmer(mean.mm.smi ~ mean.days.forelimb + (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "high density",], REML = FALSE))
lmm.test <- lmer(mean.mm.smi ~ treatment*mean.days.forelimb + (1|clutchtank), data = morph.data.mm, REML = FALSE)

test = ggeffects::ggpredict(lmm.test, terms = c("mean.days.forelimb [all]", "treatment"), type = "random", interval = "confidence")

ggplot(data=morph.data.mm, aes(x=mean.days.forelimb, y=smi, group = clutchtank, color = treatment)) +
  geom_boxplot()
  #facet_grid(cols=vars(treatment)) +
  #stat_smooth(method = "lm")


# RESPONSE TO REVIEWERS: Create Dataset from NBF Samples ----------------------------

# create dataset from individuals stored in NBF
morph.data.mm.NBF = read.csv("Database_Samples - Formalin to Ethanol Metamorphs.csv", header = TRUE, skip = 0, na.strings = "NA")

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
morph.data.mm.NBF$treatment[morph.data.mm.NBF$treatment == "control"] = "low density"
morph.data.mm.NBF$treatment = factor(morph.data.mm.NBF$treatment, levels = c("low density", "high density"))

#subset to only include RS and non-froglets and non-overflow
morph.data.mm.NBF = morph.data.mm.NBF[morph.data.mm.NBF$treatment != "overflow",]
morph.data.mm.NBF = morph.data.mm.NBF[morph.data.mm.NBF$gs == "Rana sylvatica",]
morph.data.mm.NBF = morph.data.mm.NBF[morph.data.mm.NBF$life.stage == "tail resorbed" | morph.data.mm.NBF$life.stage == "forelimb(s) emerged",]

#create unique.id.indiv
morph.data.mm.NBF$unique.id.indiv = paste(morph.data.mm.NBF$gs.code, morph.data.mm.NBF$clutch, morph.data.mm.NBF$tank.id, morph.data.mm.NBF$animal.id.num, sep = "_")

# fill developmental data column - days.forelimb for each individual
morph.data.mm.NBF$days.forelimb = NA
temp = unique(morph.data.mm.NBF$unique.id.indiv[is.na(morph.data.mm.NBF$animal.id) == FALSE])

for(i in 1:length(temp)){
  morph.data.mm.NBF$days.forelimb[morph.data.mm.NBF$unique.id.indiv == temp[i]] <- devo.data$days.forelimb[devo.data$unique.id.indiv ==  temp[i]]
}
rm(temp)

#remove NAs
morph.data.mm.NBF = morph.data.mm.NBF[is.na(morph.data.mm.NBF$sample.num)==FALSE,]
morph.data.mm.NBF = morph.data.mm.NBF[is.na(morph.data.mm.NBF$treatment)==FALSE,]

#some of these NBF individuals were part of first six but didn't survive to metamorphosis, so need to add these to first.six graph
View(morph.data.mm.NBF[morph.data.mm.NBF$animal.id <= 6,])


# ALL INDIVIDUALS WITH TAIL RESORBED - create new dataframe that combines morph.data.mm with morph.data.mm.NBF
colnames(morph.data.mm.NBF)[9] = "larv.tank.id"
morph.data.mm.NBF$clutchtank = paste(morph.data.mm.NBF$clutch, morph.data.mm.NBF$larv.tank.id, sep = "_")

morph.data.mm.ALL = rbind(
  morph.data.mm[,c(5:11, 17, 21, 35, 41)],
  morph.data.mm.NBF[,c(3:12,24)]
)

morph.data.mm.ALL = morph.data.mm.ALL[is.na(morph.data.mm.ALL$days.forelimb) == FALSE,]
morph.data.mm.ALL.TAIL = morph.data.mm.ALL[morph.data.mm.ALL$life.stage != "forelimb(s) emerged",]


#correct body length to account for fixation (equation from Shu et al. 2017. PeerJ)
morph.data.mm.NBF$svl.mm.corrected = (1.046*morph.data.mm.NBF$svl.mm) + 0.934

morph.data.mm.ALL$svl.mm.corrected[morph.data.mm.ALL$life.stage == "metamorph"] = morph.data.mm.ALL$svl.mm[morph.data.mm.ALL$life.stage == "metamorph"]

morph.data.mm.ALL.TAIL$svl.mm.corrected[morph.data.mm.ALL.TAIL$life.stage == "metamorph"] = morph.data.mm.ALL.TAIL$svl.mm[morph.data.mm.ALL.TAIL$life.stage == "metamorph"]

morph.data.mm.ALL$svl.mm.corrected[morph.data.mm.ALL$life.stage == "tail resorbed"] = (1.0946*morph.data.mm.ALL$svl.mm[morph.data.mm.ALL$life.stage == "tail resorbed"]) + 0.934

morph.data.mm.ALL.TAIL$svl.mm.corrected[morph.data.mm.ALL.TAIL$life.stage == "tail resorbed"] = (1.0946*morph.data.mm.ALL.TAIL$svl.mm[morph.data.mm.ALL.TAIL$life.stage == "tail resorbed"]) + 0.934

# RESPONSE TO REVIEWERS - PLOT DATASETS: Effect of rearing density on morphometrics AT metamorphosis considering all individuals from each tank -----------------------
first.six.mm.svl = ggplot() + 
  
  #plot first six
  geom_point(data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id.num <= 6 & morph.data.mm.ALL$life.stage == "metamorph",],
             aes(y=svl.mm, x = treatment, fill = treatment, colour = treatment),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), 
             size = 2, alpha = 0.7, show.legend = FALSE) +
  
  # geom_point(data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id <= 6& morph.data.mm.ALL$life.stage == "forelimb(s) emerged",],
  #            colour = "red",
  #            aes(y = svl.mm, x = treatment, fill = treatment),
  #            position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1),
  #            size = 2, alpha = 0.7, show.legend = FALSE) + 
  
  geom_boxplot(data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id.num <= 6 & morph.data.mm.ALL$life.stage == "metamorph",],
               aes(y=svl.mm, x = treatment, fill = treatment, colour = treatment),
               alpha = 0.6, size = 1, show.legend = FALSE, outlier.colour = "gray1", outlier.size = 1) +
  
  scale_colour_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  #add descriptive text
  #annotate(geom = "text", 
   #        label = "earliest six individuals from each tank",
    #       x=0.5, y = 15.9, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "snout-vent length (mm)", limits = c(9.5,16)) + #, breaks = seq(0.08,0.28,0.05)) +
  scale_x_discrete(name = "", na.translate = FALSE)


all.mm.svl.tail = ggplot() + 
  
  geom_point(data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id <= 6 & morph.data.mm.ALL$life.stage == "metamorph",],
             aes(y=svl.mm, x = treatment, fill = treatment, color = treatment),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), 
             size = 2, alpha = 0.7, show.legend = FALSE) +
  
  geom_point(data = morph.data.mm.NBF[morph.data.mm.NBF$animal.id > 6 & morph.data.mm.NBF$life.stage == "tail resorbed",],
             colour = "goldenrod",
             aes(y=svl.mm, x = treatment, fill = treatment),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), 
             size = 2, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "tail resorbed",],
               aes(y=svl.mm, x = treatment, fill = treatment, color = treatment),
               alpha = 0.6, size = 1, show.legend = FALSE, outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  #add descriptive text
  # annotate(geom = "text", 
  #          label = "earliest six individuals from each tank\n+all remaining individuals from each tank (gold)\nrange of sampling depth = 25% (low density)",
  #          x=0.5, y = 15.5, color = "black", size = 5, hjust=0) +
  # annotate(geom = "text", 
  #          label = "15% (high density)",
  #          x=1.44, y = 14.75, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "snout-vent length (mm)", limits = c(9.5,16)) + #, breaks = seq(0.08,0.28,0.05)) +
  scale_x_discrete(name = "", na.translate = FALSE, drop = TRUE)


all.mm.svl.tail.corrected = ggplot() + 
  
  geom_point(data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id <= 6 & morph.data.mm.ALL$life.stage == "metamorph",],
             aes(y=svl.mm.corrected, x = treatment, fill = treatment, color = treatment),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), 
             size = 2, alpha = 0.7, show.legend = FALSE) +
  
  geom_point(data = morph.data.mm.NBF[morph.data.mm.NBF$animal.id > 6 & morph.data.mm.NBF$life.stage == "tail resorbed",],
             colour = "goldenrod",
             aes(y=svl.mm.corrected, x = treatment, fill = treatment),
             position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), 
             size = 2, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "tail resorbed",],
               aes(y=svl.mm.corrected, x = treatment, fill = treatment, color = treatment),
               alpha = 0.6, size = 1, show.legend = FALSE, outlier.colour = "gray1", outlier.size = 1) +
  
  scale_color_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  scale_fill_manual(values=c(natparks.pals("Arches")[4], natparks.pals("Arches")[1])) +
  
  #add descriptive text
  annotate(geom = "text", 
           label = "earliest six individuals from each tank\n+all remaining individuals from each tank (corrected for fixation) (gold)\nrange of sampling depth = 25% (low density)",
           x=0.5, y = 15.5, color = "black", size = 5, hjust=0) +
  annotate(geom = "text", 
           label = "15% (high density)",
           x=1.44, y = 14.75, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  scale_y_continuous(name = "svl (mm)", limits = c(9.5,16)) + #, breaks = seq(0.08,0.28,0.05)) +
  scale_x_discrete(name = "", na.translate = FALSE, drop = TRUE)



ggarrange(first.six.mm.svl, all.mm.svl.tail,
          nrow = 1, ncol = 2,
          labels = c("a", "b"),
          common.legend = TRUE,
          font.label = list(size = 20, color = "black"))



# RESPONSE TO REVIEWERS - PLOT DATASETS: Effect of larval duration on morphometrics AT metamorphosis considering all individuals from each tank -----------------------

first.six.days.mm.svl = ggplot() + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, seed = 0), 
             pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",],
             colour = "black", stroke = 1.25,
             aes(y=svl.mm, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  
  geom_ribbon(data = predicted.df.morph.mm,
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.morph.mm, size = 1, color = "black", aes(x = x, y = predicted)) +
  
  scale_fill_manual(values=c("gray65", "gray65")) +

  
  # add in text
  # annotate(geom = "text",
  #          label = "earlist six individuals from each tank",
  #          x=44, y = 15.9, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic svl (mm)", limits = c(9.5,16)) +
  scale_x_continuous(name = "larval duration (days)", limits = c(43.5, 73.5), breaks = seq(44, 75, by = 5))

all.days.mm = ggplot() + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0, seed = 0), 
             pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",],
             colour = "black", stroke = 1.25,
             aes(y=svl.mm, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  
  #add in all the others
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0, seed = 0), 
             pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm.ALL[morph.data.mm.ALL$animal.id > 6 & morph.data.mm.ALL$life.stage == "tail resorbed",],
             colour = "goldenrod", stroke = 1.25,
             aes(y=svl.mm, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  
  scale_fill_manual(values=c("gray65", "gray65", "gray65")) +
  
  #add in predicted lines
  geom_ribbon(data = predicted.df.morph.mm,
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.morph.mm, size = 1, color = "black", aes(x = x, y = predicted)) +
  
  geom_ribbon(data = ALL.predicted.df.morph.mm,
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = ALL.predicted.df.morph.mm, size = 1, color = "goldenrod", aes(x = x, y = predicted)) +
  
  # # add in text
  # annotate(geom = "text",
  #          label = "earliest six individuals from each tank\n+all remaining individuals from each tank (corrected for fixation) (gold)\nrange of sampling depth = 25% (low density)",
  #          x=44, y = 15.9, color = "black", size = 5, hjust=0) +
  # annotate(geom = "text", 
  #          label = "15% (high density)",
  #          x=49, y = 15.2, color = "black", size = 5, hjust=0) +
           
  theme_bw() +
  theme(legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic svl (mm)", limits = c(9.5,16)) +
  scale_x_continuous(name = "larval duration (days)", limits = c(43.5, 73.5), breaks = seq(44, 80, by = 5))


ggarrange(first.six.days.mm.svl, all.days.mm,
          nrow = 1, ncol = 2)



# RESPONSE TO REVIEWERS - ANALYZE DATA: Effect of rearing density and larval duration on snout-vent length AT metamorphosis for FIRST SIX vs ALL---------------------

#USING ALL dataframe but only original data
# model definition 
lmm.full <- lmer(svl.mm ~ treatment*days.forelimb + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",], na.action = na.omit, REML = FALSE)

lmm.full.slopes <- lmer(svl.mm ~ treatment*days.forelimb +  (treatment||clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",], na.action = na.omit, REML = FALSE) # singular fit

lmm.nointxn <- lmer(svl.mm ~ treatment + days.forelimb + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",], na.action = na.omit, REML = FALSE)

lmm.null <- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph",], na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.full.slopes, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
final.mod = lmm.nointxn
Anova(final.mod, type = "II")
summary(final.mod)
cov2cor(vcov(final.mod)) #assess correlation matrix between fixed effects

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
leveneTest(morph.data.mm$svl.mm, morph.data.mm$treatment)
testOutliers(simulationOutput)

#post-hoc pairwise comparisons
emmeans(final.mod, pairwise ~ treatment)

# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.morph.mm <- ggpredict(final.mod, terms = c("days.forelimb [all]"), type = "random", interval = "confidence")


#using ALL dataframe and tail resorbed data
# model definition 
lmm.full <- lmer(svl.mm ~ treatment*days.forelimb + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph" | morph.data.mm.ALL$life.stage == "tail resorbed",], na.action = na.omit, REML = FALSE)

lmm.full.slopes <- lmer(svl.mm ~ treatment*days.forelimb +  (treatment||clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph" | morph.data.mm.ALL$life.stage == "tail resorbed",], na.action = na.omit, REML = FALSE) # singular fit

lmm.nointxn <- lmer(svl.mm ~ treatment + days.forelimb + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph" | morph.data.mm.ALL$life.stage == "tail resorbed",], na.action = na.omit, REML = FALSE)

lmm.null <- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.mm.ALL[morph.data.mm.ALL$life.stage == "metamorph" | morph.data.mm.ALL$life.stage == "tail resorbed",], na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.full.slopes, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
ALL.final.mod = lmm.nointxn
Anova(ALL.final.mod, type = "II")
summary(ALL.final.mod)
cov2cor(vcov(ALL.final.mod)) #assess correlation matrix between fixed effects


# check assumptions
simulationOutput <- simulateResiduals(fittedModel = ALL.final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
leveneTest(morph.data.mm$svl.mm, morph.data.mm$treatment)
testOutliers(simulationOutput)

#post-hoc pairwise comparisons
emmeans(ALL.final.mod, pairwise ~ treatment)

# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
ALL.predicted.df.morph.mm <- ggpredict(ALL.final.mod, terms = c("days.forelimb [all]"), type = "random", interval = "confidence")


# RESPONSE TO REVIEWERS - POWER ANALYSIS ---------------------

pwr.f2.test(u = 49-1, v = 397-49, f2 = 0.35, sig.level = 0.05, power = ) #49 clutch tanks, 397 observations
#numerator degrees of freedom is number of groups - 1
#denominator degrees of freedom is number of observations - number of groups

# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON MASS FOR LOW DENSITY -----------------

#define candidate models
lmm.full <- lmer(mass.g ~ days.forelimb*water.level.reduc + (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], na.action = na.omit, REML = FALSE)

lmm.nointx <- lmer(mass.g ~ days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], na.action = na.omit, REML = FALSE)

lmm.null<- lmer(mass.g ~ (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc (removing models that failed to converge or had singular fit) for non-transformed response variable
arrange(AICc(lmm.full,
             lmm.nointx, lmm.null), AICc)


# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T)
testDispersion(lmm.nointx)
testZeroInflation(lmm.nointx)
testCategorical(lmm.nointx, catPred = morph.data.mm$water.level.reduc[morph.data.mm$treatment == "low density"]) 

final.mod = lmm.nointx

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)

#create prediction dataframe to add to plots
predicted.df.morph.mm.water <- ggpredict(final.mod, terms = c("days.forelimb [all]", "water.level.reduc"), type = "random", interval = "confidence")

#plot results for metamorphic mass categorically by water level reduction
mass.waterlevel <- ggplot(data = morph.data.mm[morph.data.mm$treatment == "low density",], aes(y=mass.g, x = water.level.reduc, fill = water.level.reduc, color = water.level.reduc)) + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0, seed = 1), size = 1.5, alpha = 0.7, show.legend = FALSE) +
  
  geom_boxplot(alpha = 0.6, size = 1, show.legend = FALSE, aes(fill = water.level.reduc), outlier.colour = "gray1", outlier.size = 1, color = "black") +
  
  stat_summary(fun = mean, geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.1, size = 1, alpha=1, color = "black") +
  
  stat_summary(fun=mean, geom="point", color = "black", pch=21, size=6, stroke = 1, show.legend = FALSE) +
  
  scale_fill_manual(values=c("gray10", "goldenrod")) +
  scale_color_manual(values=c("gray10", "goldenrod")) +
  
  theme_bw() +
  theme(#legend.position = c(0.02, 0.98),
        #legend.justification = c("left", "top"),
        legend.title = element_blank(),
        #legend.text = element_text(size=10),
        axis.text.x = element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  geom_signif(comparisons = list(c("no", "yes")), 
              map_signif_level=TRUE,
              annotations = "ns",
              size = 1,
              tip_length = 0,
              vjust = -0.3,
              textsize = 6,
              color = "black") + 
  
  scale_y_continuous(name = "metamorphic mass (g)", limits = c(0.07,0.29), breaks = seq(0.08,0.28,0.05)) +
  scale_x_discrete(name = "", labels = c("not reduced", "reduced"))

#plot results for relationship between metamorphic mass and larval duration based on water level reduction
mass.larvdur.waterlevel <- ggplot() + 
  
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0, seed = 0), 
             pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm[morph.data.mm$water.level.reduc == "no" & morph.data.mm$treatment == "low density",],
             colour = "black", stroke = 1.25,
             aes(y=mass.g, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  
  #add in all the others
  geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0, dodge.width = 0, seed = 0), 
             pch = 21, alpha = 0.7, size = 3,
             data = morph.data.mm[morph.data.mm$water.level.reduc == "yes" & morph.data.mm$treatment == "low density",],
             colour = "goldenrod", stroke = 1.25,
             aes(y=mass.g, x = days.forelimb, fill = treatment), show.legend = FALSE) +
  
  scale_fill_manual(values=c("gray65", "gray65", "gray65")) +
  
  #add in predicted lines
  geom_ribbon(data = predicted.df.morph.mm.water[predicted.df.morph.mm.water$group == "no",],
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.morph.mm.water[predicted.df.morph.mm.water$group == "no",], size = 1, color = "black", aes(x = x, y = predicted)) +
  
  geom_ribbon(data = predicted.df.morph.mm.water[predicted.df.morph.mm.water$group == "yes",],
              alpha = 0.4, 
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.morph.mm.water[predicted.df.morph.mm.water$group == "yes",], size = 1, color = "goldenrod", aes(x = x, y = predicted)) +
  
  # # add in text
  # annotate(geom = "text",
  #          label = "earliest six individuals from each tank\n+all remaining individuals from each tank (corrected for fixation) (gold)\nrange of sampling depth = 25% (low density)",
  #          x=44, y = 15.9, color = "black", size = 5, hjust=0) +
  # annotate(geom = "text", 
  #          label = "15% (high density)",
  #          x=49, y = 15.2, color = "black", size = 5, hjust=0) +
  
  theme_bw() +
  theme(legend.background=element_rect(fill = alpha("white", 0)), #make legend background transparent
        axis.text.x=element_text(size=20, color = "black"), 
        axis.text.y=element_text(size=20, color = "black"), 
        axis.title.x=element_text(size=20, color = "black"), 
        axis.title.y = element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "metamorphic mass (g)", limits = c(0.07,0.29), breaks = seq(0.08,0.28,0.05)) +
  scale_x_continuous(name = "larval duration (days)", limits = c(43.5, 73.5), breaks = seq(44, 80, by = 5))

ggarrange(mass.waterlevel, mass.larvdur.waterlevel,
          labels = c("a","b"),
          font.label = list(size = 20, color = "black")
  )


# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON DAYS FORELIMB FOR LOW DENSITY -----------------

#define candidate models
lmm.nointx <- lmer(days.forelimb ~ water.level.reduc + (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], na.action = na.omit, REML = FALSE)

lmm.null<- lmer(days.forelimb ~ (1|clutchtank), data = morph.data.mm[morph.data.mm$treatment == "low density",], na.action = na.omit, REML = FALSE)

# non-nested model selection using AICc (removing models that failed to converge or had singular fit) for non-transformed response variable
arrange(AICc(lmm.nointx, lmm.null), AICc)

# check assumptions of best-fit model with non-transformed response variable
simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T)
testDispersion(lmm.nointx)
testZeroInflation(lmm.nointx)
testCategorical(lmm.nointx, catPred = morph.data.mm$water.level.reduc[morph.data.mm$treatment == "low density"]) 

final.mod = lmm.nointx

# estimates from best-supported model
Anova(final.mod, type = "II")
summary(final.mod)


# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON METAMORPHIC STATUS FOR LOW DENSITY, CALCULATED FROM ALL SEEDED INDIVIDUALS -----------------
#define candidate models
glmm.full.01 <- glmer(status ~ week*water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longform[survi.data.tadexp.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.full.01.slopes <- glmer(status ~ week*water.level.reduc + (week||clutchtank), data=survi.data.tadexp.longform[survi.data.tadexp.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ week + water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longform[survi.data.tadexp.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01.slopes <- glmer(status ~ week + water.level.reduc + (week||clutchtank), data=survi.data.tadexp.longform[survi.data.tadexp.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.01 <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tadexp.longform[survi.data.tadexp.longform$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc (removing models that failed to converge or had singular fit)
arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
             glmm.nointxn.01, glmm.nointxn.01.slopes,
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(
  rownames(arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
                        glmm.nointxn.01, 
                        glmm.null.01), AICc))[1])))

#check assumptions
simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tadexp.longform$week[survi.data.tadexp.longform$treatment == "low density"])
plotResiduals(glmm.full.01, form = survi.data.tadexp.longform$water.level.reduc[survi.data.tadexp.longform$treatment == "low density"])
testDispersion(glmm.full.01)


# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant


# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON METAMORPHIC STATUS FOR LOW DENSITY, CALCULATED FROM ALL SEEDED INDIVIDUALS BUT NOT USING LONGFORM -----------------
#define candidate models
glmm.full.01 <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ week*water.level.reduc + (1|clutchtank), data=survi.data.tadexp[survi.data.tadexp$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), weights = (seed.num-leth.samp.num.cumul))

glmm.nointxn.01 <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ week + water.level.reduc + (1|clutchtank), data=survi.data.tadexp[survi.data.tadexp$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), weights = (seed.num-leth.samp.num.cumul))

glmm.null.01 <- glmer(metamorph.num.cumul/(seed.num-leth.samp.num.cumul) ~ 1 + (1|clutchtank), data=survi.data.tadexp[survi.data.tadexp$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), weights = (seed.num-leth.samp.num.cumul))

# non-nested model selection using AICc (removing models that failed to converge or had singular fit)
arrange(AICc(glmm.full.01,
             glmm.nointxn.01,
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(
  rownames(arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
                        glmm.nointxn.01, 
                        glmm.null.01), AICc))[1])))

#check assumptions
simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tadexp$week[survi.data.tadexp$treatment == "low density"])
plotResiduals(glmm.full.01, form = survi.data.tadexp$water.level.reduc[survi.data.tadexp$treatment == "low density"])
testDispersion(glmm.full.01)


# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant

# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON METAMORPHIC STATUS FOR LOW DENSITY, CALCULATED FROM SURVIVING INDIVIDUALS -----------------
#define candidate models
glmm.full.01 <- glmer(status ~ week*water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longformsurvi[survi.data.tadexp.longformsurvi$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ week + water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longformsurvi[survi.data.tadexp.longformsurvi$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.01 <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tadexp.longformsurvi[survi.data.tadexp.longformsurvi$treatment == "low density",], na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc (removing models that failed to converge or had singular fit)
arrange(AICc(glmm.full.01, 
             glmm.nointxn.01, 
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(
  rownames(arrange(AICc(glmm.full.01, glmm.full.01.slopes, 
                        glmm.nointxn.01, 
                        glmm.null.01), AICc))[1])))

#check assumptions
simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tadexp.longformsurvi$week[survi.data.tadexp.longformsurvi$treatment == "low density"])
plotResiduals(glmm.full.01, form = survi.data.tadexp.longformsurvi$water.level.reduc[survi.data.tadexp.longformsurvi$treatment == "low density"])
testDispersion(glmm.full.01)


# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant



# RESPONSE TO REVIEWERS - EFFECT OF WATER LEVEL ON METAMORPHIC STATUS FOR BOTH DENSITIES, CALCULATED FROM ALL SEEDED INDIVIDUALS -----------------
#define candidate models
glmm.full.01 <- glmer(status ~ week + treatment*water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn.01 <- glmer(status ~ week + treatment + water.level.reduc + (1|clutchtank), data=survi.data.tadexp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null.01 <- glmer(status ~ 1 + (1|clutchtank), data=survi.data.tadexp.longform, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# non-nested model selection using AICc (removing models that failed to converge or had singular fit)
arrange(AICc(glmm.full.01, 
             glmm.nointxn.01, 
             glmm.null.01), AICc)

final.mod = eval(parse(text = paste(
  rownames(arrange(AICc(glmm.full.01, 
                        glmm.nointxn.01, 
                        glmm.null.01), AICc))[1])))

#check assumptions
simulateResiduals(fittedModel = glmm.full.01, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(glmm.full.01, form = survi.data.tadexp.longform$treatment)
plotResiduals(glmm.full.01, form = survi.data.tadexp.longform$week)
plotResiduals(glmm.full.01, form = survi.data.tadexp.longform$water.level.reduc)
testDispersion(glmm.full.01)


# estimates from best-supported model
summary(final.mod)
Anova(final.mod, type = "III", contrasts=list(topic=contr.sum, sys=contr.sum)) #type 3 because interaction term significant