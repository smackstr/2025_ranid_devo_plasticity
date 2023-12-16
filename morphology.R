# this script creates figures of the morphology results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on morphology metrics, 2. effects of larval duration on morphology metrics) across three life stages (1. larval period, 2. at metamorphosis, 3. up to 3 months post-metamorphosis)

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

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for metamorph morphometrics, juvenile morphometrics, and developmental timing data
morph.data.tad = read.csv("Database_Morphometrics - Tadpole Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv = read.csv("Database_Morphometrics - Froglet_Toadlet Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

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

# subset devo dataset to only include first six individuals and non-overflow individuals
devo.data = devo.data[devo.data$first.six == "yes",]
devo.data = devo.data[devo.data$treatment != "overflow",]

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
morph.data.juv$post.mm.weeks = factor(morph.data.juv$post.mm.weeks, ordered = TRUE, levels = c("1-2", "4-6", "5-7", "8-10", "11-12", "12-14", "16-18")) #set weeks as ordered factor

morph.data.tad$tmw.cm = as.numeric(morph.data.tad$tmw.cm)

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
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")

# need to change juvenile tanks with "extra" in their name to be just the name of the tank
morph.data.juv$unique.id = paste(morph.data.juv$gs.code, morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")
devo.data$unique.id.juv = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")

# create unique individual id for devo.data and morph.data.mm, since these are the two we have individual-level data for
morph.data.mm$unique.id.indiv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, morph.data.mm$animal.id, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, devo.data$animal.id, sep = "_")

#remove any NAs that have occurred
devo.data = devo.data[is.na(devo.data$gs) == FALSE,]
morph.data.tad = morph.data.tad[is.na(morph.data.tad$gs) == FALSE,]


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


# create column so can plot post.mm.weeks on numeric scale
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "1-2"] = mean(c(1,2))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "4-6"] = mean(c(4,6))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "8-10"] = mean(c(8,10))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "12-14"] = mean(c(12,14))
morph.data.juv$post.mm.weeks.num[morph.data.juv$post.mm.weeks == "16-18"] = mean(c(16,18))


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


# create column to store categorical developmental data (i.e. early, mid, late developers)
morph.data.mm$devo.cat = NA
morph.data.mm$devo.cat[morph.data.mm$days.forelimb <= (mean(devo.data$days.forelimb, na.rm = TRUE) - sd(devo.data$days.forelimb, na.rm = TRUE))] = "early"
morph.data.mm$devo.cat[morph.data.mm$days.forelimb >= (mean(devo.data$days.forelimb, na.rm = TRUE) + sd(devo.data$days.forelimb, na.rm = TRUE))] = "late"   
morph.data.mm$devo.cat[is.na(morph.data.mm$devo.cat) == TRUE] = "mid"


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.mm.juv -----------------------
# create column to store the developmental data
morph.data.mm.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.mm.juv$unique.id))){
  morph.data.mm.juv$mean.days.forelimb[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.mm.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical developmental data (i.e. early, mid, late developers)
morph.data.mm.juv$devo.cat = NA

for(i in 1:length(unique(morph.data.mm.juv$unique.id))){
  
  if(morph.data.mm.juv$mean.days.forelimb[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]][1] <= mean(devo.data$days.forelimb, na.rm = TRUE) - (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    morph.data.mm.juv$devo.cat[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]] = "early"
  }
  
  if(morph.data.mm.juv$mean.days.forelimb[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]][1] >= mean(devo.data$days.forelimb, na.rm = TRUE) + (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    morph.data.mm.juv$devo.cat[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i]] = "late"
  }
}
morph.data.mm.juv$devo.cat[is.na(morph.data.mm.juv$devo.cat) == TRUE] = "mid"


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.juv -----------------------
# create column to store the developmental data
morph.data.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id.juv ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical developmental data (i.e. early, mid, late developers)
morph.data.juv$devo.cat = NA

for(i in 1:length(unique(morph.data.juv$unique.id))){
  
  if(morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]][1] <= mean(devo.data$days.forelimb, na.rm = TRUE) - (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    morph.data.juv$devo.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] = "early"
  }
  
  if(morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]][1] >= mean(devo.data$days.forelimb, na.rm = TRUE) + (sd(devo.data$days.forelimb, na.rm = TRUE)/2)){
    morph.data.juv$devo.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] = "late"
  }
}
morph.data.juv$devo.cat[is.na(morph.data.juv$devo.cat) == TRUE] = "mid"


# COMPILE DATASETS: Add growth rate column (growth.mass.g & growth.svl.mm) to morph.data.tad -----------------------
morph.data.tad = morph.data.tad[is.na(morph.data.tad$tl.cm) == FALSE,]

morph.data.tad$unique.id.date = paste(morph.data.tad$unique.id, morph.data.tad$week, sep = "_")
morph.data.tad$growth.tl.cm = NA
morph.data.tad$growth.bl.cm = NA
morph.data.tad$date.measured =  as.Date(morph.data.tad$date.photo, format = "%Y-%m-%d")


growth.data.tad = morph.data.tad %>%
  group_by(week, gs, gs.code, clutch, treatment, treatment.code, larv.tank.id) %>%
  filter(is.na(larv.tank.id) == FALSE) %>%
  filter(week > 2) %>%
  summarise(n = n())


growth.data.tad$unique.id = paste(growth.data.tad$gs.code, growth.data.tad$clutch, growth.data.tad$larv.tank.id, sep = "_")
growth.data.tad$unique.id.date = paste(growth.data.tad$gs.code, growth.data.tad$clutch, growth.data.tad$larv.tank.id, growth.data.tad$week, sep = "_")
growth.data.tad$growth.tl.cm = NA
growth.data.tad$growth.bl.cm = NA

for(i in 1:length(unique(morph.data.tad$unique.id))){
  temp1 = morph.data.tad[morph.data.tad$unique.id == unique(morph.data.tad$unique.id)[i],]
  for (j in 2:length(unique(temp1$unique.id.date))){
    temp2 = temp1[temp1$unique.id.date == unique(temp1$unique.id.date)[j],] #current week
    temp3 = temp1[temp1$unique.id.date == unique(temp1$unique.id.date)[j-1],] #previous week
    
    growth.data.tad$growth.tl.cm[growth.data.tad$unique.id.date == temp2$unique.id.date[1]] = (mean(temp2$tl.cm, na.rm = TRUE) - mean(temp3$tl.cm, na.rm = TRUE)) / as.numeric((max(temp2$date.measured, na.rm = TRUE) - min(temp3$date.measured, na.rm = TRUE))[1])
    
    growth.data.tad$growth.bl.cm[growth.data.tad$unique.id.date == temp2$unique.id.date[1]] = (mean(temp2$bl.cm, na.rm = TRUE) - mean(temp3$bl.cm, na.rm = TRUE)) / as.numeric((max(temp2$date.measured, na.rm = TRUE) - min(temp3$date.measured, na.rm = TRUE))[1])
  }
}


# COMPILE DATASETS: Create growth rate database for juvenile growth -----------------------
# morph.data.mm.juv$unique.id.date = paste(morph.data.mm.juv$unique.id, morph.data.mm.juv$post.mm.weeks, sep = "_")
# morph.data.mm.juv$growth.mass.g = NA
# morph.data.mm.juv$date.measured =  as.Date(morph.data.mm.juv$date.measured, format = "%Y-%m-%d")

growth.data.mm.juv = morph.data.mm.juv %>%
  group_by(post.mm.weeks, gs, gs.code, clutch, treatment, treatment.code, juv.tank.id, mean.days.forelimb, devo.cat) %>%
  filter(is.na(juv.tank.id) == FALSE) %>%
  filter(post.mm.weeks != "0") %>%
  summarise(n = n())

growth.data.mm.juv$unique.id = paste(growth.data.mm.juv$gs.code, growth.data.mm.juv$clutch, growth.data.mm.juv$juv.tank.id, sep = "_")
growth.data.mm.juv$unique.id.date = paste(growth.data.mm.juv$gs.code, growth.data.mm.juv$clutch, growth.data.mm.juv$juv.tank.id, growth.data.mm.juv$post.mm.weeks, sep = "_")
growth.data.mm.juv$growth.mass.g = NA
growth.data.mm.juv$growth.svl.mm = NA


# fill mean.mm.mass.g for growth.data.mm.juv -- currently calculating mean based on first six mass
growth.data.mm.juv$mean.mm.mass.g = NA
for(i in 1:length(unique(growth.data.mm.juv$unique.id))){
  growth.data.mm.juv$mean.mm.mass.g[growth.data.mm.juv$unique.id == unique(growth.data.mm.juv$unique.id)[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(growth.data.mm.juv$unique.id)[i]], na.rm = TRUE)
}

# create column so can plot post.mm.weeks on numeric scale
growth.data.mm.juv$post.mm.weeks.num = NA
growth.data.mm.juv$post.mm.weeks.num[growth.data.mm.juv$post.mm.weeks == "1-2"] = mean(c(1,2))
growth.data.mm.juv$post.mm.weeks.num[growth.data.mm.juv$post.mm.weeks == "4-6"] = mean(c(4,6))
growth.data.mm.juv$post.mm.weeks.num[growth.data.mm.juv$post.mm.weeks == "8-10"] = mean(c(8,10))
growth.data.mm.juv$post.mm.weeks.num[growth.data.mm.juv$post.mm.weeks == "12-14"] = mean(c(12,14))
growth.data.mm.juv$post.mm.weeks.num[growth.data.mm.juv$post.mm.weeks == "16-18"] = mean(c(16,18))


#create unique id for morph.data.mm.juv$unique.id that has sampling weeks included
morph.data.mm.juv$unique.id.date = paste(morph.data.mm.juv$unique.id, morph.data.mm.juv$post.mm.weeks, sep = "_")

# fill other columns with growth
for(i in 1:length(unique(morph.data.mm.juv$unique.id))){
  temp1 = morph.data.mm.juv[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id)[i],]
  for (j in 2:length(unique(temp1$unique.id.date))){
    temp2 = temp1[temp1$unique.id.date == unique(temp1$unique.id.date)[j],] #current week
    temp3 = temp1[temp1$unique.id.date == unique(temp1$unique.id.date)[j-1],] #previous week
    
    #set date to be a date
    temp2$date.measured = as.Date(temp2$date.measured)
    temp3$date.measured = as.Date(temp3$date.measured)
    
    growth.data.mm.juv$growth.mass.g[growth.data.mm.juv$unique.id.date == temp2$unique.id.date[1]] = (mean(temp2$mass.g, na.rm = TRUE) - mean(temp3$mass.g, na.rm = TRUE)) / as.numeric(max(temp2$date.measured, na.rm = TRUE) - min(temp3$date.measured, na.rm = TRUE))
    growth.data.mm.juv$growth.svl.mm[growth.data.mm.juv$unique.id.date == temp2$unique.id.date[1]] = (mean(temp2$svl.mm, na.rm = TRUE) - mean(temp3$svl.mm, na.rm = TRUE)) / as.numeric(max(temp2$date.measured, na.rm = TRUE) - min(temp3$date.measured, na.rm = TRUE))
  }
}


# create column to store mean.mm.smi, which represents average smi at metamorphosis for each juvenile tank
growth.data.mm.juv$mean.mm.smi = NA

# fill mean.mm.smi for growth.data.mm.juv
for(i in 1:length(unique(growth.data.mm.juv$unique.id))){
  growth.data.mm.juv$mean.mm.smi[growth.data.mm.juv$unique.id == unique(growth.data.mm.juv$unique.id)[i]] <- mean(morph.data.mm$smi[morph.data.mm$unique.id ==  unique(growth.data.mm.juv$unique.id)[i]], na.rm = TRUE)
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

# create column to store categorical size data (i.e. small, medium, large tanks)
morph.data.mm.juv$mass.cat = NA

for(i in 1:length(unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE]))){
  
  if(morph.data.mm.juv$mean.mm.mass.g[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]][1] <= mean(morph.data.mm$mass.g, na.rm = TRUE) - (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    morph.data.mm.juv$mass.cat[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]] = "small"
  }
  
  if(morph.data.mm.juv$mean.mm.mass.g[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]][1] >= mean(morph.data.mm$mass.g, na.rm = TRUE) + (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    morph.data.mm.juv$mass.cat[morph.data.mm.juv$unique.id == unique(morph.data.mm.juv$unique.id[is.na(morph.data.mm.juv$juv.tank.id) == FALSE])[i]] = "large"
  }
}
morph.data.mm.juv$mass.cat[is.na(morph.data.mm.juv$mass.cat) == TRUE & is.na(morph.data.mm.juv$juv.tank.id) == FALSE] = "med"


# COMPILE DATASETS: Add Metamorphosis size (mean mass at mm) to morph.data.juv -----------------------

# create column to store the developmental data, mean mass at mm represents average mass at metamorphosis for each juvenile tank
morph.data.juv$mean.mm.mass.g = NA

# fill mean.mm.mass.g for morph.data.juv -- currently calculating mean based on first six mass
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.mm.mass.g[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(morph.data.mm$mass.g[morph.data.mm$unique.id ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}

# create column to store categorical size data (i.e. small, medium, large tanks)
morph.data.juv$mass.cat = NA

for(i in 1:length(unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE]))){
  
  if(morph.data.juv$mean.mm.mass.g[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]][1] <= mean(morph.data.mm$mass.g, na.rm = TRUE) - (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    morph.data.juv$mass.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]] = "small"
  }
  
  if(morph.data.juv$mean.mm.mass.g[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]][1] >= mean(morph.data.mm$mass.g, na.rm = TRUE) + (sd(morph.data.mm$mass.g, na.rm = TRUE)/2)){
    morph.data.juv$mass.cat[morph.data.juv$unique.id == unique(morph.data.juv$unique.id[is.na(morph.data.juv$juv.tank.id) == FALSE])[i]] = "large"
  }
}
morph.data.juv$mass.cat[is.na(morph.data.juv$mass.cat) == TRUE & is.na(morph.data.juv$juv.tank.id) == FALSE] = "med"
morph.data.juv$mass.cat = factor(morph.data.juv$mass.cat, levels = c("small", "med", "large"))


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


#COMPILE DATASETS: Create clutchtank for later analysis   -----------------------
morph.data.tad$clutchtank = paste(morph.data.tad$clutch, morph.data.tad$larv.tank.id, sep="_")
morph.data.mm$clutchtank = paste(morph.data.mm$clutch, morph.data.mm$larv.tank.id, sep="_")
morph.data.mm.juv$clutchtank = paste(morph.data.mm.juv$clutch, morph.data.mm.juv$juv.tank.id, sep="_")
morph.data.juv$clutchtank = paste(morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep="_")
growth.data.tad$clutchtank = paste(growth.data.tad$clutch, growth.data.tad$larv.tank.id, sep="_")



# COMPILE DATASETS: Create metrics vector to feed into plots and analyses later  -----------------------
metrics.tad = colnames(morph.data.tad)[17:22]
metrics.mm = colnames(morph.data.mm)[17:21]
metrics.mm.juv = colnames(morph.data.mm.juv)[c(11:15,21)]
metrics.juv = colnames(morph.data.juv)[c(24:28,33)]



# ANALYZE DATA: Effect of rearing density and age on larval morphometrics TOTAL LENGTH ---------------------

# model definition - using glmer with log-link function to keep original units.
lmm.full <- lmer(tl.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log <- lmer(log(tl.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log <- lmer(log(tl.cm)~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly <- lmer(log(tl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(tl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx <- lmer(tl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(tl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointx, lmm.null, test="Chisq")
anova(lmm.full, lmm.full.log, lmm.full.log.log, lmm.full.log.poly, lmm.full.log.poly.slopes)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.tad$treatment) 
testCategorical(simulationOutput, catPred = morph.data.tad$week) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)


#examine the plots for for each tank, treatment, and clutch, regress tl.cm by week using linear and then using linear and you can see how much better the glm is and how it does look like distinct slopes for high density vs low density
ggplot(morph.data.tad, 
       aes(x = week, y = tl.cm, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.tad, method = "lm", aes(x=week, y=tl.cm), inherit.aes = F, se = F, color="black") #linear

ggplot(morph.data.tad, 
       aes(x = log(week), y = tl.cm, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.tad, aes(x=log(week), y=tl.cm), inherit.aes = F, se = F, color="black") #quadratic

# Final model
final.mod = lmm.full.log.poly
Anova(final.mod, type = "II")
summary(final.mod)

predicted.df.tad.tl <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density and age on larval morphometrics BODY LENGTH ---------------------

# model definition - using glmer with log-link function to keep original units.
lmm.full <- lmer(bl.cm ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log <- lmer(log(bl.cm) ~ treatment*week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.log <- lmer(log(bl.cm)~ treatment*log(week) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly <- lmer(log(bl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (1|clutchtank), data = morph.data.tad, na.action = na.omit)
lmm.full.log.poly.slopes <- lmer(log(bl.cm)~ treatment*poly(week, degree = 2, raw = FALSE) + (week||clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.nointx <- lmer(bl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

lmm.null <- lmer(bl.cm ~ treatment + week + (1|clutchtank), data = morph.data.tad, na.action = na.omit)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointx, lmm.null, test="Chisq")
anova(lmm.full, lmm.full.log, lmm.full.log.log, lmm.full.log.poly, lmm.full.log.poly.slopes)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full.log.poly, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.tad$treatment) 
testCategorical(simulationOutput, catPred = morph.data.tad$week) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)


#examine the plots for for each tank, treatment, and clutch, regress tl.cm by week using linear and then using linear and you can see how much better the glm is and how it does look like distinct slopes for high density vs low density
ggplot(morph.data.tad, 
       aes(x = week, y = bl.cm, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.tad, method = "lm", aes(x=week, y=bl.cm), inherit.aes = F, se = F, color="black") #linear

ggplot(morph.data.tad, 
       aes(x = log(week), y = bl.cm, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.tad, aes(x=log(week), y=bl.cm), inherit.aes = F, se = F, color="black") #quadratic

# Final model
final.mod = lmm.full.log.poly
Anova(final.mod, type = "II")
summary(final.mod)

predicted.df.tad.bl <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")



# ANALYZE DATA: Effect of rearing density and age on larval growth TOTAL LENGTH ----------
growth.data.tad = growth.data.tad[is.na(growth.data.tad$growth.tl.cm) == FALSE,]
#model definition - using glmer with log-link function to keep original units.
lmm.full <- lmer(growth.tl.cm ~ treatment*week + (1|clutchtank), data = growth.data.tad, REML = FALSE)

lmm.nointx <- lmer(growth.tl.cm ~ treatment + week + (1|clutchtank), data = growth.data.tad, REML = FALSE)

lmm.null <- lmer(growth.tl.cm ~ (1|clutchtank), data = growth.data.tad, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointx, lmm.null, test="Chisq")

# Final Model
final.mod = lmm.full
Anova(final.mod, type = "II")
summary(final.mod)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = growth.data.tad$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)

predicted.df.tad.growth.tl <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density and age on larval growth BODY LENGTH ----------
growth.data.tad = growth.data.tad[is.na(growth.data.tad$growth.bl.cm) == FALSE,]
#model definition - using glmer with log-link function to keep original units.
lmm.full <- lmer(growth.bl.cm ~ treatment*week + (1|clutchtank), data = growth.data.tad, REML = FALSE)
lmm.full.poly <- lmer(growth.bl.cm ~ treatment*poly(week, degree = 2, raw=FALSE) + (1|clutchtank), data = growth.data.tad, REML = FALSE)
lmm.full.poly.slopes <- lmer(growth.bl.cm ~ treatment*poly(week, degree = 2, raw=FALSE) + (week||clutchtank), data = growth.data.tad, REML = FALSE)
lmm.full.slopes <- lmer(growth.bl.cm ~ treatment*week + (week||clutchtank), data = growth.data.tad, REML = FALSE)

lmm.nointx <- lmer(growth.bl.cm ~ treatment + week + (1|clutchtank), data = growth.data.tad, REML = FALSE)

lmm.null <- lmer(growth.bl.cm ~ (1|clutchtank), data = growth.data.tad, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointx, lmm.null, test="Chisq")
anova(lmm.full, lmm.full.poly, lmm.full.poly.slopes, lmm.null)

# Final Model
final.mod = lmm.full.poly
Anova(final.mod, type = "II")
summary(final.mod)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full.poly, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = growth.data.tad$treatment) 
testCategorical(simulationOutput, catPred = growth.data.tad$week) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)

predicted.df.tad.growth.bl <- ggeffects::ggpredict(final.mod, terms = c("week[all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on laterality of forelimb emergence FIRST SIX ---------------------
# model definition - 

glmm.full <- glmer(details.forelimb.binary ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.nointxn <- glmer(details.forelimb.binary ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmm.null <- glmer(details.forelimb.binary ~ (1|clutchtank), data = morph.data.mm, na.action = na.omit, family = binomial, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

# model selection using likelihood ratio test
anova(glmm.full, glmm.nointxn, glmm.null, test="Chisq")

#models did not differ from null
final.mod = glmm.null

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = glmm.null, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on mass AT metamorphosis for FIRST SIX---------------------

# model definition 

lmm.full <- lmer(mass.g ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.full.log <- lmer(log(mass.g) ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(mass.g ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn.log <- lmer(log(mass.g) ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(log(mass.g) ~ (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")
anova(lmm.nointxn)

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
testOutliers(simulationOutput)
leveneTest(morph.data.mm$mass.g~morph.data.mm$treatment)

# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.morph.mm <- ggpredict(final.mod, terms = c("days.forelimb [all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on snout-vent length AT metamorphosis for FIRST SIX---------------------
# model definition 
lmm.full <- lmer(svl.mm ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(svl.mm ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.log <- lmer(log(svl.mm) ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
final.mod = lmm.nointxn
Anova(lmm.nointxn, type = "II")
summary(lmm.nointxn)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
leveneTest(morph.data.mm$svl.mm, morph.data.mm$treatment)
testOutliers(simulationOutput)

#post-hoc pairwise comparisons
emmeans(final.mod, pairwise ~ treatment)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on forearm length AT metamorphosis for FIRST SIX ---------------------
# model definition 
lmm.full <- lmer(r.forelimb.mm ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(r.forelimb.mm ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(r.forelimb.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
final.mod = lmm.nointxn
Anova(final.mod, type = "II")
summary(final.mod)
cov2cor(vcov(lmm.nointxn)) #assess correlation matrix between fixed effects

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
leveneTest(morph.data.mm$r.forelimb.mm, morph.data.mm$treatment)
testOutliers(simulationOutput)

#based on DHARMaoutputs we have an issue with heteroskedasticity across treatment. A rule of thumb is that linear models are fairly robust to heterogeneity of variance so long as the maximum variance is no more than 4× greater than the minimum variance, which is demonstrated by the following ratio. So if this ratio is <4, you are likely okay with running the linear model
var(morph.data.mm$r.forelimb.mm[morph.data.mm$treatment == "high density"], na.rm=T) / var(morph.data.mm$r.forelimb.mm[morph.data.mm$treatment == "low density"], na.rm=T)

#but we can also check using a model that accounts for heteroskedasticity and allows the error variance to depend on fixed effects in addition to random effects to vary using nlme package; uses a generalized least squares regression (gls), which is a form of weighted regression
lme.nointxn <- lme(fixed = r.forelimb.mm ~ treatment + scale(days.forelimb) + water.level.reduc + mass.g,
                random = (~1|clutchtank),
                data = morph.data.mm,
                na.action = na.omit)

lme.nointxn.var <- lme(fixed = r.forelimb.mm ~ treatment + scale(days.forelimb) + water.level.reduc + mass.g,
                    random = (~1|clutchtank),
                    data = morph.data.mm,
                    weights = varIdent(form=~1|treatment),
                    na.action = na.omit)

# check assumptions
par(mfrow = c(1,1))
plot(lme.nointxn, resid(., type = "p") ~ fitted(.) | treatment, abline = 0 )
plot(lme.nointxn.var, resid(., type = "p") ~ fitted(.) | treatment, abline = 0 )

# model selection using likelihood ratio test to determine if accounting for heterskedasticity significantly improves gls model fit
anova(lme.nointxn, lme.nointxn.var)

# Final Model: 
# check if different results if white.adjust set to true to allow for heteroskedasticity
car::Anova(lmm.nointxn, white.adjust=TRUE, type = "II")
car::Anova(lmm.nointxn, white.adjust=FALSE, type = "II")

final.mod = lme.nointxn.var
Anova(final.mod, type = "II")
summary(final.mod)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on tibia length AT metamorphosis for FIRST SIX ---------------------
# model definition 

lmm.full <- lmer(r.tibia.mm ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(r.tibia.mm ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(r.tibia.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

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
testOutliers(simulationOutput)
leveneTest(morph.data.mm$r.tibia.mm, morph.data.mm$treatment)

#post-hoc pairwise comparisons
emmeans(final.mod, pairwise ~ treatment)
emmeans(final.mod, pairwise ~ water.level.reduc)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on thigh length AT metamorphosis for FIRST SIX ---------------------
# model definition 
lmm.full <- lmer(r.thigh.mm ~ treatment*scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(r.thigh.mm ~ treatment + scale(days.forelimb) + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(r.thigh.mm ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

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
testOutliers(simulationOutput)
leveneTest(morph.data.mm$r.thigh.mm, morph.data.mm$treatment)

#post-hoc pairwise comparisons
emmeans(final.mod, pairwise ~ treatment)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on scaled mass index (SMI = body condition) AT metamorphosis for FIRST SIX ---------------------
# model definition - using glmer with log-link function to keep mass in original units. create new column to store because ggpredict doesn't work to recalculate formulaic response variables

#optimize random effects with most saturated model

lmm.full.1 <- lmer(smi ~ treatment*days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = TRUE)
lmm.full.1.slopes <- lmer(smi ~ treatment*days.forelimb + water.level.reduc + (treatment|clutchtank), data = morph.data.mm, na.action = na.omit, REML = TRUE)

AICc(lmm.full.1,lmm.full.1.slopes)

#optimize fixed effects with optimal random effects structure
lmm.full.1 <- lmer(smi ~ treatment*days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.full.1.log <- lmer(log(smi) ~ treatment*days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(smi ~ treatment + days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nointxn.log <- lmer(log(smi) ~ treatment + days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nodays <- lmer(smi ~ days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.nodays.log <- lmer(log(smi) ~ days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(smi ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)
lmm.null.log <- lmer(log(smi) ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using log-likelihood - "Fit no more parameters than is necessary. If two or more models fit the data almost equally well, prefer the simpler one
anova(lmm.full.1, lmm.nointxn, lmm.nodays, lmm.null)
anova(lmm.full.1.log, lmm.nointxn.log, lmm.nodays.log, lmm.null.log)

# Final Model
final.mod = lmm.null

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full.1.log, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)
leveneTest(morph.data.mm$smi, group = morph.data.mm$treatment)


# ANALYZE DATA: Effect of rearing density, water level reduction, and larval duration on scaled mass index (SMI = body condition) AT metamorphosis for FIRST SIX BUT NOW USING TANK MEANS ---------------------
# model definition - using glmer with log-link function to keep mass in original units. create new column to store because ggpredict doesn't work to recalculate formulaic response variables

#optimize random effects with most saturated model

lmm.full.1 <- lmer(mean.mm.smi ~ treatment*mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = TRUE)
lmm.full.1.slopes <- lmer(mean.mm.smi ~ treatment*mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = TRUE)

AICc(lmm.full.1,lmm.full.1.slopes)

#optimize fixed effects with optimal random effects structure
lmm.full.1 <- lmer(mean.mm.smi ~ treatment*mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.full.1.log <- lmer(log(mean.mm.smi) ~ treatment*mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn <- lmer(mean.mm.smi ~ treatment + mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nointxn.log <- lmer(log(mean.mm.smi) ~ treatment + mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nodays <- lmer(mean.mm.smi ~ mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.nodays.log <- lmer(log(mean.mm.smi) ~ mean.days.forelimb + water.level.reduc + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null <- lmer(mean.mm.smi ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

lmm.null.log <- lmer(log(mean.mm.smi) ~ 1 + (1|clutchtank), data = morph.data.mm, na.action = na.omit, REML = FALSE)

# model selection using log-likelihood - "Fit no more parameters than is necessary. If two or more models fit the data almost equally well, prefer the simpler one
anova(lmm.full.1, lmm.nointxn, lmm.nodays, lmm.null)
anova(lmm.full.1.log, lmm.nointxn.log, lmm.nodays.log, lmm.null.log)

# Final Model
final.mod = lmm.full.1.log
Anova(final.mod, type = "II")
summary(final.mod)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full.1.log, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.mm$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)
leveneTest(morph.data.mm$mean.mm.smi, group = morph.data.mm$treatment)


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on mass after metamorphosis ---------------------
# model definition
# QUESTION: should I use juv.tank.id nested within clutch if mean.mm.mass.g is synonymous with juv.tank.id because it's not changing
lmm.full <- lmer(mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(mass.g ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(mass.g ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(mass.g ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly3 <- lmer(mass.g ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(mass.g ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null.log<- lmer(log(mass.g) ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested models
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)


# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))

cov2cor(vcov(final.mod)) #assess correlation matrix between fixed effects

emmeans(final.mod, pairwise ~ treatment)

cor.test(morph.data.juv$mean.mm.mass.g, morph.data.juv$mean.days.forelimb)


# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.morph.juv <- ggpredict(final.mod, terms = c("post.mm.weeks.num [all]"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on mass after metamorphosis BUT NOW USING TANK-LEVEL MEANS ---------------------
# data examination
ggplot(morph.data.juv.summ, 
       aes(x = post.mm.weeks.num, y = mean.mass.g, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.juv.summ, method = "lm", aes(x=post.mm.weeks.num, y=mean.mass.g), inherit.aes = F, se = F, color="black") #linear

ggplot(morph.data.juv.summ, 
       aes(x = log(post.mm.weeks.num), y = mean.mass.g, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.juv.summ, aes(x=log(post.mm.weeks.num), y=mean.mass.g), inherit.aes = F, se = F, color="black") #log

morph.data.juv.summ$poly.post.mm.weeks.num = morph.data.juv.summ$post.mm.weeks.num^2
ggplot(morph.data.juv.summ, 
       aes(x = post.mm.weeks.num^2, y = mean.mass.g, group = clutchtank, color = clutchtank)) +
  #geom_line() +
  stat_summary(fun.y=mean, geom="line", size = 0.5, aes(color = clutchtank, group = clutchtank)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=2, aes(fill=clutchtank)) +
  facet_wrap(~treatment) +
  geom_smooth(data = morph.data.juv.summ, aes(x=post.mm.weeks.num^2, y = mean.mass.g, inherit.aes = F, se = F, color="black")) #poly


# model definition
# QUESTION: should I use juv.tank.id nested within clutch if mean.mm.mass.g is synonymous with juv.tank.id because it's not changing
lmm.full <- lmer(mean.mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(mean.mass.g ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(mean.mass.g ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(mean.mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(mean.mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(mean.mass.g) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(mean.mass.g ~ (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.null.log <- lmer(log(mean.mass.g) ~ (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested model structure
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv.summ$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv.summ$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)

# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))

cov2cor(vcov(final.mod)) #assess correlation matrix between fixed effects

emmeans(final.mod, pairwise ~ treatment)

cor.test(morph.data.juv$mean.mm.mass.g, morph.data.juv$mean.days.forelimb)


# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.morph.juv.summ <- ggpredict(final.mod, terms = c("post.mm.weeks.num [all]"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on SVL after metamorphosis ---------------------
# model definition
lmm.full <- lmer(svl.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(svl.mm ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(svl.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(svl.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly3 <- lmer(svl.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(svl.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(svl.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(svl.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(svl.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null.log<- lmer(log(svl.mm) ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested models
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)


# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on FOREARM after metamorphosis ---------------------
# model definition
lmm.full <- lmer(r.forelimb.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.forelimb.mm ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.forelimb.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.forelimb.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly3 <- lmer(r.forelimb.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(r.forelimb.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.forelimb.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(r.forelimb.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.forelimb.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null.log<- lmer(log(r.forelimb.mm) ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested models
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)


# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on TIBIA LENGTH after metamorphosis ---------------------
# model definition
lmm.full <- lmer(r.tibia.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.tibia.mm ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.tibia.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.tibia.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly3 <- lmer(r.tibia.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(r.tibia.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.tibia.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(r.tibia.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.tibia.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null.log<- lmer(log(r.tibia.mm) ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested models
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)


# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on THIGH LENGTH after metamorphosis ---------------------
# model definition
lmm.full <- lmer(r.thigh.mm ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.1 <- lmer(r.thigh.mm ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(r.thigh.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly2 <- lmer(r.thigh.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.poly3 <- lmer(r.thigh.mm ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log <- lmer(log(r.thigh.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly2 <- lmer(log(r.thigh.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 2) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn.log.poly3 <- lmer(log(r.thigh.mm) ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + poly(post.mm.weeks.num, degree = 3) + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(r.thigh.mm ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null.log<- lmer(log(r.thigh.mm) ~ (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn.1, lmm.nointxn, lmm.null)

# model selection using AICc because non-nested models
AICc(lmm.nointxn, lmm.nointxn.poly2, lmm.nointxn.poly3, lmm.null)
AICc(lmm.nointxn.log, lmm.nointxn.log.poly2, lmm.nointxn.log.poly3, lmm.null.log)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointxn.log.poly3, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$treatment)
plotResiduals(lmm.nointxn.log.poly3, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)


# Final Model
final.mod = lmm.nointxn.log.poly3
Anova(final.mod, type = "II")
summary(final.mod)
exp(fixef(final.mod))


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on scaled mass index (body condition) after metamorphosis ---------------------

# model definition
lmm.full <- lmer(smi ~ treatment*scale(mean.days.forelimb)*mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.intxn1 <- lmer(smi ~ treatment + scale(mean.days.forelimb)*mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(smi ~ treatment + scale(mean.days.forelimb) + mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(smi ~ 1 + (1|clutchtank), data = morph.data.juv, na.action = na.omit, REML=FALSE)


# model selection using likelihood ratio test
anova(lmm.full, lmm.intxn1, lmm.nointxn, lmm.null, test="Chisq")

#model selection using AICc since non-nested models
AICc(lmm.full, lmm.intxn1, lmm.nointxn, lmm.treat, lmm.null)

# Final Model
final.mod = lmm.full
Anova(final.mod, type = "II")
summary(final.mod)

#three-way interaction is significant so helpful to examine relationship separately for hd and ld
lmm.full.ld <- lmer(smi ~ scale(mean.days.forelimb)*mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv[morph.data.juv$treatment == "low density",], na.action = na.omit, REML=FALSE)
lmm.full.hd <- lmer(smi ~ scale(mean.days.forelimb)*mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv[morph.data.juv$treatment == "high density",], na.action = na.omit, REML=FALSE)

Anova(lmm.full.ld, type = "II")
Anova(lmm.full.hd, type = "II")

summary(lmm.full.ld)
summary(lmm.full.hd)

vcov(lmm.full) |> cov2cor() #assess correlations between fixed effects
emmeans(lmm.full, pairwise ~ treatment)

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
plotResiduals(lmm.full, form = morph.data.juv$treatment)
plotResiduals(lmm.full, form = morph.data.juv$post.mm.weeks.num)
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.juv$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)

#post-hoc pairwise comparisons
emmeans(lmm.full, pairwise ~ treatment)

# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.bodycond.juv <- ggpredict(final.mod, terms = c("post.mm.weeks.num [all]", "treatment"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of rearing density, metamorphic mass, and larval duration on scaled mass index (body condition) after metamorphosis BUT NOW USING TANK-LEVEL MEANS ---------------------

# model definition
lmm.full <- lmer(mean.smi ~ treatment*scale(mean.days.forelimb)*mean.mm.smi+ post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.intxn1 <- lmer(mean.smi ~ treatment + scale(mean.days.forelimb)*mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.nointxn <- lmer(mean.smi ~ treatment + scale(mean.days.forelimb) + mean.mm.smi + post.mm.weeks.num + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(mean.smi ~ 1 + (1|clutchtank), data = morph.data.juv.summ, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.intxn1, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
final.mod = lmm.nointxn
Anova(final.mod, type = "II")
summary(final.mod)

vcov(final.mod) |> cov2cor() #assess correlations between fixed effects

# check assumptions
simulationOutput <- simulateResiduals(fittedModel = final.mod, quantreg=T, plot = T) #provides summary of model fitting tests

testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.juv.summ$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)

#post-hoc pairwise comparisons
emmeans(final.mod, pairwise ~ treatment)

# create dataframe of predicted values that can be plotted on ggplot later
# predictions using ggeffects (suggested by Susan Durham)
predicted.df.bodycond.juv.summ <- ggpredict(final.mod, terms = c("post.mm.weeks.num [all]"), type = "random", interval = "confidence")


# ANALYZE DATA: Effect of treatment and larval duration on mean smi to show that larval duration influences mean smi ---------------------
# model definition

lmm.full <- lmer(mean.mm.smi ~ treatment*scale(days.forelimb) + (1|clutch), data = morph.data.mm, na.action = na.omit)

lmm.nointxn <- lmer(log(mean.mm.smi) ~ treatment + days.forelimb + (1|clutch), data = morph.data.mm, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(log(mean.mm.smi) ~ (1|clutch), data = morph.data.mm, na.action = na.omit, REML=FALSE)


# model selection using likelihood ratio test
anova(lmm.full, lmm.nointxn, lmm.null, test="Chisq")

# Final Model
Anova(lmm.full, type = "III")
summary(lmm.full)
# R's default is to fit a series of polynomial functions or contrasts to the levels of the variable. Variables with .L, .Q, and .C and ^4 are, respectively, the coefficients for the ordered factor coded with linear, quadratic, cubic, and quadratic contrasts. The first is linear (.L), the second is quadratic (.Q), the third is cubic (.C), and so on. R will fit one fewer polynomial functions than the number of available levels. The interpretation of a particular beta test is then generalized to: Which contrasts contribute significantly to explain any differences between levels in your dependent variable? Because the weeks.L predictor is significant and negative, this suggests a linear decreasing trend in logit across weeks, and because the Year.Q predictor is significant and negative, this suggests a deacceleration trend is detectable in the pattern of logits across years. Third order polynomials model jerk (rate of change of an object's acceleration over time), and fourth order polynomials model jounce (a.k.a., snap). However, I would stop interpreting around this order and higher because it quickly becomes nonsensical to practical folk.

#back-transform fixed effects
fixef(lmm.full)

emmeans(lmm.full, pairwise ~ treatment)


# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.full, quantreg=T, plot = T) #provides summary of model fitting tests
testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testOutliers(simulationOutput)



# ANALYZE DATA: Effect of rearing density and larval duration on mass after metamorphosis BUT NOW USING CATEGORY FOR METAMORPHIC MASS ---------------------

# model definition - using glmer with log-link function to keep mass in original units.
lmm.full <- lmer(log(mass.g) ~ treatment*scale(mean.days.forelimb)*post.mm.weeks + mass.cat + (1|clutch), data = morph.data.juv[is.na(morph.data.juv$juv.tank.id)==FALSE,], na.action = na.omit, REML=FALSE)

# create dataframe of predicted values that can be plotted on ggplot later. predictions using ggeffects (suggested by Susan Durham)
predicted.df.juv <- ggpredict(lmm.full, terms = c("post.mm.weeks", "mass.cat"), type = "random")
# create column so can plot x on numeric scale
predicted.df.juv$x.num[predicted.df.juv$x == "1-2"] = mean(c(1,2))
predicted.df.juv$x.num[predicted.df.juv$x == "4-6"] = mean(c(4,6))
predicted.df.juv$x.num[predicted.df.juv$x == "8-10"] = mean(c(8,10))
predicted.df.juv$x.num[predicted.df.juv$x == "12-14"] = mean(c(12,14))
predicted.df.juv$x.num[predicted.df.juv$x == "16-18"] = mean(c(16,18))


# ANALYZE DATA: Effect of rearing density and larval duration on scaled mass index after metamorphosis BUT NOW USING CATEGORY FOR SMI ---------------------

lmm.full <- lmer(smi ~ treatment*devo.cat*smi.cat + post.mm.weeks + (1|clutch), data = morph.data.juv[is.na(morph.data.juv$juv.tank.id)==FALSE,], na.action = na.omit, REML=FALSE)

#predict for later ggplot
predicted.df.bodycond.juv <- ggpredict(lmm.full, terms = c("post.mm.weeks", "smi.cat"), type = "random", interval = "confidence")
# create column so can plot x on numeric scale
predicted.df.bodycond.juv$x.num[predicted.df.bodycond.juv$x == "1-2"] = mean(c(1,2))
predicted.df.bodycond.juv$x.num[predicted.df.bodycond.juv$x == "4-6"] = mean(c(4,6))
predicted.df.bodycond.juv$x.num[predicted.df.bodycond.juv$x == "8-10"] = mean(c(8,10))
predicted.df.bodycond.juv$x.num[predicted.df.bodycond.juv$x == "12-14"] = mean(c(12,14))
predicted.df.bodycond.juv$x.num[predicted.df.bodycond.juv$x == "16-18"] = mean(c(16,18))



# ANALYZE DATA: Effect of rearing density and larval duration on growth rate after metamorphosis ---------------------
# model definition - using glmer with log-link function to keep mass in original units.
growth.data.mm.juv$post.mm.weeks = factor(growth.data.mm.juv$post.mm.weeks, ordered = TRUE)

# model definition
lmm.full <- lmer(growth.mass.g ~ treatment*scale(mean.days.forelimb)*mean.mm.mass.g + mean.mm.smi + post.mm.weeks + (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

lmm.intxn1 <- lmer(growth.mass.g  ~ treatment*scale(mean.days.forelimb) + mean.mm.mass.g + mean.mm.smi + post.mm.weeks + (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

lmm.intxn2 <- lmer(growth.mass.g  ~ treatment + scale(mean.days.forelimb)*mean.mm.mass.g + mean.mm.smi + post.mm.weeks + (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

lmm.intxn3 <- lmer(growth.mass.g  ~ treatment*mean.mm.mass.g + scale(mean.days.forelimb) + mean.mm.smi + post.mm.weeks + (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

lmm.nointx <- lmer(growth.mass.g ~ treatment + scale(mean.days.forelimb) + mean.mm.mass.g + mean.mm.smi + post.mm.weeks + (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

lmm.null<- lmer(growth.mass.g  ~ (1|clutch:juv.tank.id), data = growth.data.mm.juv, na.action = na.omit, REML=FALSE)

# model selection using likelihood ratio test
anova(lmm.full, lmm.intxn1, lmm.intxn2, lmm.intxn3, lmm.nointx, lmm.null, test="Chisq")

# Final Model
Anova(lmm.nointx, type = "II")
summary(lmm.nointx)
# R's default is to fit a series of polynomial functions or contrasts to the levels of the variable. Variables with .L, .Q, and .C and ^4 are, respectively, the coefficients for the ordered factor coded with linear, quadratic, cubic, and quadratic contrasts. The first is linear (.L), the second is quadratic (.Q), the third is cubic (.C), and so on. R will fit one fewer polynomial functions than the number of available levels. The interpretation of a particular beta test is then generalized to: Which contrasts contribute significantly to explain any differences between levels in your dependent variable? Because the weeks.L predictor is significant and negative, this suggests a linear decreasing trend in logit across weeks, and because the Year.Q predictor is significant and negative, this suggests a deacceleration trend is detectable in the pattern of logits across years. Third order polynomials model jerk (rate of change of an object's acceleration over time), and fourth order polynomials model jounce (a.k.a., snap). However, I would stop interpreting around this order and higher because it quickly becomes nonsensical to practical folk.
exp(fixef(lmm.nointx)[-1] + fixef(lmm.nointx)[1])


# check assumptions
simulationOutput <- simulateResiduals(fittedModel = lmm.nointx, quantreg=T, plot = T) #provides summary of model fitting tests

testDispersion(simulationOutput) #tests for over- and under-dispersion
testZeroInflation(simulationOutput) #tests if more zeroes than expected
testCategorical(simulationOutput, catPred = morph.data.juv$treatment) #tests residuals against a categorical predictor to assess homogeneity of variance; heteroscedasticity means that there is a systematic dependency of the dispersion / variance on another variable in the model...it means that the level of over/underdispersion depends on another parameter.
testOutliers(simulationOutput)

#predict for later ggplot
predicted.df.growth.juv <- ggpredict(lmm.nointx, terms = c("post.mm.weeks"), type = "random", interval = "confidence")
# create column so can plot x on numeric scale
predicted.df.growth.juv$x.num[predicted.df.growth.juv$x == "1-2"] = mean(c(1,2))
predicted.df.growth.juv$x.num[predicted.df.growth.juv$x == "4-6"] = mean(c(4,6))
predicted.df.growth.juv$x.num[predicted.df.growth.juv$x == "8-10"] = mean(c(8,10))
predicted.df.growth.juv$x.num[predicted.df.growth.juv$x == "12-14"] = mean(c(12,14))
predicted.df.growth.juv$x.num[predicted.df.growth.juv$x == "16-18"] = mean(c(16,18))


# ANALYZE DATA: Effect of rearing density, treatment, and larval duration on snout-vent length at and after metamorphosis ---------------------

# model definition - REML set to FALSE because comparing models with different fixed effects using hypothesis test (likelihood ratio test). REML assumes that all fixed effects are the same in comparison models but ML does not. setting random effect as juvenile tank id nested within clutch because want to control for correlation within those groups but not necessarily interested in defining the effect of being in those groups. have to specify data = na.omit(morph.data.mm.juv) because some NAs in mass.g so model fitted without mass.g will have issues when we compare to other models
glmer.full <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb)*post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

glmer.nointx <- glmer(svl.mm ~ treatment + scale(mean.days.forelimb) + post.mm.weeks + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))

glmer.null<- glmer(svl.mm ~ (1|post.mm.weeks) + (1|clutch:juv.tank.id), data = na.omit(morph.data.mm.juv), na.action = na.omit, family = gaussian(link = 'log'))


# model selection using likelihood ratio test
anova(glmer.full, glmer.nointx, glmer.null, test="Chisq")

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


# ANALYZE DATA: Correlation between morphology metrics during larval period -----------------------

# create a dataframe of the correlation between each metric and mass for each rearing density and clutch and post-metamorphic age
for(i in 1:length(unique(morph.data.tad$week))){
  for(j in 1:length(unique(morph.data.tad$treatment))){
    temp1 = morph.data.tad[morph.data.tad$week == unique(morph.data.tad$week)[i] & morph.data.tad$treatment == unique(morph.data.tad$treatment)[j],]
    
    corr.bl <- cor.test(temp1$tl.cm, temp1$bl.cm, 
                         method = "pearson", na.action = na.omit)
    corr.tail <- cor.test(temp1$tl.cm, temp1$tail.cm, 
                              method = "pearson", na.action = na.omit)
    corr.hw <- cor.test(temp1$tl.cm, temp1$hw.cm, 
                           method = "pearson", na.action = na.omit)
    corr.bw <- cor.test(temp1$tl.cm, temp1$bw.cm, 
                           method = "pearson", na.action = na.omit)
    corr.tmw <- cor.test(temp1$tl.cm, temp1$tmw.cm, 
                           method = "pearson", na.action = na.omit)
    
    temp2 = data.frame(metric = c("bl.cm", "tail.cm", "hw.cm", "bw.cm", "tmw.cm"),
                       treatment = unique(temp1$treatment)[1],
                       week = unique(temp1$week)[1],
                       type = c("bl.cm", "tail.cm", "hw.cm", "bw.cm", "tmw.cm"),
                       t = c(corr.bl$statistic, corr.tail$statistic, corr.hw$statistic, corr.bw$statistic, corr.tmw$statistic),
                       df = c(corr.bl$parameter, corr.tail$parameter, corr.hw$parameter, corr.bw$parameter, corr.tmw$parameter),
                       p.value = c(corr.bl$p.value, corr.tail$p.value, corr.hw$p.value, corr.bw$p.value, corr.tmw$p.value),
                       estimate = c(corr.bl$estimate, corr.tail$estimate, corr.hw$estimate, corr.bw$estimate, corr.tmw$estimate),
                       conf.int2.5 = c(corr.bl$conf.int[1], corr.tail$conf.int[1], corr.hw$conf.int[1], corr.bw$conf.int[1], corr.tmw$conf.int[1]),
                       conf.int97.5 = c(corr.bl$conf.int[2], corr.tail$conf.int[2], corr.hw$conf.int[2], corr.bw$conf.int[2], corr.tmw$conf.int[2])
                       )
    
    if(i == 1 & j == 1){
      morph.data.corr.tad = temp2}else{
        morph.data.corr.tad = rbind(morph.data.corr.tad, temp2)
      }
    rm(temp1, temp2)
  }
}


# ANALYZE DATA: Correlation between morphology metrics at and after metamorphosis -----------------------

# create a dataframe of the correlation between each metric and mass for each rearing density and clutch and post-metamorphic age
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
                       post.mm.weeks = unique(temp1$post.mm.weeks)[1],
                       post.mm.sampling = unique(temp1$post.mm.sampling)[1],
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
glm.full <- glm(estimate ~ treatment*post.mm.weeks + type , data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

glm.nointxn <- glm(estimate ~ treatment+post.mm.weeks + type , data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

glm.notreat <- glm(estimate ~ post.mm.weeks + type , data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

glm.notype <- glm(estimate ~ treatment + post.mm.weeks, data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

glm.noweeks <- glm(estimate ~ treatment + type, data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

glm.null <- glm(estimate ~ 1, data = morph.data.corr, na.action = na.omit, family = gaussian(link = 'log'))

# model selection using likelihood ratio test
anova(glm.full, glm.nointxn, glm.notreat, glm.notype, glm.noweeks, glm.null, test="Chisq")
logLik(glm.full)
logLik(glm.notype) #first model significantly different from the full but logLik is maximized for the full so going with that one

# Check Model Assumptions
check_model(glm.full)

# Final Model
summary(glm.full)
Anova(glm.full, type = "III") #gives information for each factor; should use type = "II" when interaction terms are NOT significant and thus not included in the final model; should use type = "III" when interaction terms are significant and thus included in the final model
confint(glm.full)


# PLOT DATASETS: Effect of rearing density on morphometrics during larval stage -----------------------
#morph.data.tad = morph.data.tad[morph.data.tad$week.range != "",]
#morph.data.tad$week.range = factor(morph.data.tad$week.range, ordered = TRUE, levels = c("2","3","4-5","6-7","8-9","10"))

ggplot(data = morph.data.tad, aes(y=tl.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
  geom_jitter(size = 2.5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment)) +
  #geom_boxplot(alpha = 0.6, size = 0.75) +
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
  scale_y_continuous(name = "total length (cm)") +
  scale_x_discrete(name = "week")

#plot tanks as linked through time
ggplot(data = morph.data.tad, aes(y=tl.cm, x = as.numeric(week, ordered = TRUE, levels = c(2,3,4,5,6,7,8,9,10)), color = treatment)) + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  #facet_grid(rows = vars(clutch)) +
  geom_jitter(size = 2.5, alpha = 0.6, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="line", size = 0.8, aes(color = treatment, group = larv.tank.id)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = larv.tank.id)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment, group = larv.tank.id)) +
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
  scale_y_continuous(name = "total length (cm)", breaks = c(1.5,2,2.5,3,3.5,4,4.5)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

#plot tanks as linked through time with predicted values from model
plot.morph.tad.tl.ld <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_jitter(data = morph.data.tad[morph.data.tad$treatment == "low density",], width=0.4, size = 2, alpha = 0.7, pch = 21, aes(y=tl.cm, x = week, fill = treatment, color = treatment)) +
  geom_smooth(data = morph.data.tad[morph.data.tad$treatment == "low density",], 
              fun.y=mean, geom="line", size = 0.3, se=F, 
              aes(y=tl.cm, x = week, group = larv.tank.id, color = treatment)) +
  stat_summary(data = morph.data.tad[morph.data.tad$treatment == "low density",], 
               fun.y=mean, geom="point", pch=21, size=5, 
               aes(y=tl.cm, x = week, fill=treatment, color = treatment, group = larv.tank.id)) +
  geom_ribbon(data = predicted.df.tad.tl[predicted.df.tad.tl$group == "low density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.tl[predicted.df.tad.tl$group == "low density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3])) +
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
  scale_y_continuous(name = "total length (cm)", breaks = c(1.5,2,2.5,3,3.5,4,4.5)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

plot.morph.tad.tl.hd <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_jitter(data = morph.data.tad[morph.data.tad$treatment == "high density",], width=0.4, size = 2, alpha = 0.7, pch = 21, aes(y=tl.cm, x = week, fill = treatment, color = treatment)) +
  geom_smooth(data = morph.data.tad[morph.data.tad$treatment == "high density",], 
              fun.y=mean, geom="line", size = 0.3, se=F, 
              aes(y=tl.cm, x = week, group = larv.tank.id, color = treatment)) +
  stat_summary(data = morph.data.tad[morph.data.tad$treatment == "high density",], 
               fun.y=mean, geom="point", pch=21, size=5, 
               aes(y=tl.cm, x = week, fill=treatment, color = treatment, group = larv.tank.id)) +
  geom_ribbon(data = predicted.df.tad.tl[predicted.df.tad.tl$group == "high density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.tl[predicted.df.tad.tl$group == "high density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=natparks.pals("BryceCanyon")[1]) +
  scale_fill_manual(values=natparks.pals("BryceCanyon")[1]) +
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
  scale_y_continuous(name = "total length (cm)", breaks = c(1.5,2,2.5,3,3.5,4,4.5)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

ggplot(data = morph.data.tad, aes(y=bl.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "body length (cm)") +
  scale_x_discrete(name = "week")


#plot tanks as linked through time with predicted values from model
plot.morph.tad.bl.ld <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_jitter(data = morph.data.tad[morph.data.tad$treatment == "low density",], width=0.4, size = 2, alpha = 0.7, pch = 21, aes(y=bl.cm, x = week, fill = treatment, color = treatment)) +
  geom_smooth(data = morph.data.tad[morph.data.tad$treatment == "low density",], 
              fun.y=mean, geom="line", size = 0.3, se=F, 
              aes(y=bl.cm, x = week, group = larv.tank.id, color = treatment)) +
  stat_summary(data = morph.data.tad[morph.data.tad$treatment == "low density",], 
               fun.y=mean, geom="point", pch=21, size=5, 
               aes(y=bl.cm, x = week, fill=treatment, color = treatment, group = larv.tank.id)) +
  geom_ribbon(data = predicted.df.tad.bl[predicted.df.tad.bl$group == "low density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.bl[predicted.df.tad.bl$group == "low density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3])) +
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
  scale_y_continuous(name = "body length (cm)", limits = c(0.5,1.75), breaks = seq(0.5,1.75,0.25)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

plot.morph.tad.bl.hd <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_jitter(data = morph.data.tad[morph.data.tad$treatment == "high density",], width=0.4, size = 2, alpha = 0.7, pch = 21, aes(y=bl.cm, x = week, fill = treatment, color = treatment)) +
  geom_smooth(data = morph.data.tad[morph.data.tad$treatment == "high density",], 
              fun.y=mean, geom="line", size = 0.3, se=F, 
              aes(y=bl.cm, x = week, group = larv.tank.id, color = treatment)) +
  stat_summary(data = morph.data.tad[morph.data.tad$treatment == "high density",], 
               fun.y=mean, geom="point", pch=21, size=5, 
               aes(y=bl.cm, x = week, fill=treatment, color = treatment, group = larv.tank.id)) +
  geom_ribbon(data = predicted.df.tad.bl[predicted.df.tad.bl$group == "high density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.bl[predicted.df.tad.bl$group == "high density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=natparks.pals("BryceCanyon")[1]) +
  scale_fill_manual(values=natparks.pals("BryceCanyon")[1]) +
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
  scale_y_continuous(name = "body length (cm)", limits = c(0.5,1.75), breaks = seq(0.5,1.75,0.25)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))



ggplot(data = morph.data.tad, aes(y=bl.cm/tl.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "body length/total length") +
  scale_x_discrete(name = "week")

ggplot(data = morph.data.tad, aes(y=tail.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "tail length (cm)") +
  scale_x_discrete(name = "week")


ggplot(data = morph.data.tad, aes(y=tail.cm/tl.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "tail length/total length") +
  scale_x_discrete(name = "week")


ggplot(data = morph.data.tad, aes(y=hw.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "head width (cm)") +
  scale_x_discrete(name = "week")


ggplot(data = morph.data.tad, aes(y=bw.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "body width (cm)") +
  scale_x_discrete(name = "week")

ggplot(data = morph.data.tad, aes(y=bw.cm/bl.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "body width/body length (body condition)") +
  scale_x_discrete(name = "week")


ggplot(data = morph.data.tad, aes(y=tmw.cm, x = factor(week), color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "tail muscle width (cm)") +
  scale_x_discrete(name = "week")


# PLOT DATASETS: Effect of rearing density on growth rate during larval stage -----------------------
ggplot(data = growth.data.tad, aes(y=growth.tl.cm, x = week, color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
  geom_jitter(size = 2.5, alpha = 1, show.legend = FALSE) +
  stat_summary(fun.y=mean, geom="line", size = 1.2, aes(color = treatment, group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.7, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, aes(fill=treatment)) +
  #geom_boxplot(alpha = 0.6, size = 0.75) +
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
  scale_y_continuous(name = "growth in total length (cm/day)") +
  scale_x_continuous(name = "week")

#plot tanks as linked through time and density in paired plot
ggplot(data = growth.data.tad, aes(y=growth.tl.cm, x = week, color = treatment)) + 
  facet_grid(rows = vars(treatment)) +
  geom_hline(yintercept=0, size = 0.5, alpha = 1, linetype = 2, show.legend = FALSE) +
  geom_point(size = 0.9, alpha = 0.6, show.legend = FALSE, aes(color = treatment, group = larv.tank.id)) +
  geom_line(alpha = 0.6, size = 0.6, aes(color = treatment, group = larv.tank.id)) +
  stat_summary(fun.y=mean, geom="line", alpha = 0.6, size = 1, color = "black", aes(group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.6, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, alpha = 0.9, na.rm = TRUE, aes(fill=treatment, group = treatment)) +
  #geom_boxplot(alpha = 0.6, size = 0.75) +
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
  scale_y_continuous(name = "growth in total length (cm/day)") +
  scale_x_discrete(name = "week")

#plot tanks as linked through time and density in shared plot
ggplot(data = growth.data.tad, aes(y=growth.tl.cm, x = as.numeric(week, ordered = FALSE), color = treatment)) + 
  #facet_grid(rows = vars(treatment)) +
  geom_hline(yintercept=0, size = 0.5, alpha = 1, color = "black", linetype = 2, show.legend = FALSE) +
  geom_vline(xintercept=(44/49)*5, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_point(size = 0.9, alpha = 0.6, show.legend = FALSE, aes(color = treatment, group = larv.tank.id)) +
  geom_line(alpha = 0.6, size = 0.3, aes(color = treatment, group = larv.tank.id)) +
  stat_summary(fun.y=mean, geom="line", alpha = 0.6, size = 1, color = "black", aes(group = treatment)) +
  stat_summary(fun = mean,
               geom = "errorbar",
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)), #plotting +1 se
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)), #plotting -1 se
               width=0.07, size = 0.7, colour="black", alpha=0.6, aes(group = treatment)) +
  stat_summary(fun.y=mean, geom="point", color = "black", pch=21, size=5, alpha = 0.9, na.rm = TRUE, aes(fill=treatment, group = treatment)) +
  #geom_boxplot(alpha = 0.6, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x = element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "growth in total length (cm/day)") +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8), labels = unique(factor(growth.data.tad$week)))

#plot tanks as linked through time with predicted values from model
plot.morph.tad.growth.tl.ld <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_point(data = growth.data.tad[growth.data.tad$treatment == "low density",], size = 2, alpha = 0.7, pch = 21, aes(y=growth.tl.cm, x = week, fill = treatment, color = treatment, group= larv.tank.id)) +
  geom_line(data = growth.data.tad[growth.data.tad$treatment == "low density",], 
              size = 0.3,
              aes(y=growth.tl.cm, x = week, group = larv.tank.id, color = treatment)) +
  geom_ribbon(data = predicted.df.tad.growth.tl[predicted.df.tad.growth.tl$group == "low density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.growth.tl[predicted.df.tad.growth.tl$group == "low density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3])) +
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
  scale_y_continuous(name = "growth in total length (cm)", limits = c(-0.02, 0.06), breaks = seq(-0.02,0.06,0.02)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

plot.morph.tad.growth.tl.hd <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_point(data = growth.data.tad[growth.data.tad$treatment == "high density",], size = 2, alpha = 0.7, pch = 21, aes(y=growth.tl.cm, x = week, fill = treatment, color = treatment, group= larv.tank.id)) +
  geom_line(data = growth.data.tad[growth.data.tad$treatment == "high density",], 
            size = 0.3,
            aes(y=growth.tl.cm, x = week, group = larv.tank.id, color = treatment)) +
  geom_ribbon(data = predicted.df.tad.growth.tl[predicted.df.tad.growth.tl$group == "high density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.growth.tl[predicted.df.tad.growth.tl$group == "high density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[1])) +
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
  scale_y_continuous(name = "growth in total length (cm)", limits = c(-0.02, 0.06), breaks = seq(-0.02,0.06,0.02)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))


ggplot(data = growth.data.tad, aes(y=growth.bl.cm, x = week, color = treatment)) + 
  #facet_grid(rows = vars(clutch)) +
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
  scale_y_continuous(name = "growth in body length (cm/day)") +
  scale_x_discrete(name = "week")


#playing around with connecting tanks with lines
ggplot(data = growth.data.tad, aes(y=growth.tl.cm, x = factor(week), color = treatment, group = larv.tank.id)) + 
  #facet_grid(rows = vars(clutch)) +
  geom_line(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, seed = 0.2), size = 0.5, alpha = 1, show.legend = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0, seed = 0.2), size = 2.5, alpha = 1, show.legend = FALSE) +
  #geom_boxplot(alpha = 0.6, size = 0.75) +
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
  scale_y_continuous(name = "growth in total length (cm/day)") +
  scale_x_discrete(name = "week")


#plot tanks as linked through time with predicted values from model
plot.morph.tad.growth.bl.ld <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_point(data = growth.data.tad[growth.data.tad$treatment == "low density",], size = 2, alpha = 0.7, pch = 21, aes(y=growth.bl.cm, x = week, fill = treatment, color = treatment, group= larv.tank.id)) +
  geom_line(data = growth.data.tad[growth.data.tad$treatment == "low density",], 
            size = 0.3,
            aes(y=growth.bl.cm, x = week, group = larv.tank.id, color = treatment)) +
  geom_ribbon(data = predicted.df.tad.growth.bl[predicted.df.tad.growth.bl$group == "low density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.growth.bl[predicted.df.tad.growth.bl$group == "low density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3])) +
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
  scale_y_continuous(name = "growth in body length (cm)", limits = c(-0.02, 0.06), breaks = seq(-0.02,0.06,0.02)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))

plot.morph.tad.growth.bl.hd <- ggplot() + 
  geom_vline(xintercept=(44/49)*7, size = 0.5, alpha = 1, color = "black", linetype = 3, show.legend = FALSE) + #add metamorphosis start date
  geom_point(data = growth.data.tad[growth.data.tad$treatment == "high density",], size = 2, alpha = 0.7, pch = 21, aes(y=growth.bl.cm, x = week, fill = treatment, color = treatment, group= larv.tank.id)) +
  geom_line(data = growth.data.tad[growth.data.tad$treatment == "high density",], 
            size = 0.3,
            aes(y=growth.bl.cm, x = week, group = larv.tank.id, color = treatment)) +
  geom_ribbon(data = predicted.df.tad.growth.bl[predicted.df.tad.growth.bl$group == "high density",],
              alpha = 0.7, color = "black", show.legend = FALSE,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group)) +
  geom_line(data = predicted.df.tad.growth.bl[predicted.df.tad.growth.bl$group == "high density",],
            alpha = 1, size = 1.8, color = "black", show.legend = FALSE,
            mapping = aes(x = x, y = predicted, group = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[1])) +
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
  scale_y_continuous(name = "growth in body length (cm)", limits = c(-0.02, 0.06), breaks = seq(-0.02,0.06,0.02)) +
  scale_x_continuous(name = "week", breaks = c(1,2,3,4,5,6,7,8,9,10))


# PLOT DATASETS: Effect of rearing density on morphometrics AT metamorphosis for first six individuals from each tank -----------------------
plot.morph1 <- ggplot(data = morph.data.mm, aes(y=mass.g, x = clutch, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  geom_boxplot(alpha = 0.6, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
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
  scale_y_continuous(name = "mass (g) at metamorphosis") +
  scale_x_discrete(name = "clutches separated")

ggplot(data = morph.data.mm, aes(y=svl.mm, x = clutch, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  geom_boxplot(alpha = 0.6, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
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
  scale_y_continuous(name = "svl (cm) at metamorphosis") +
  scale_x_discrete(name = "clutches separated")

ggplot(data = morph.data.mm, aes(y=smi, x = clutch, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 1, show.legend = FALSE) +
  geom_boxplot(alpha = 0.6, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
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
  scale_y_continuous(name = "scaled mass index at metamorphosis") +
  scale_x_discrete(name = "clutches separated")


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

ggplot(data = morph.data.corr.tad, aes(y=estimate, x = factor(week), color = treatment)) + 
  facet_grid(rows = vars(metric)) +
  geom_point(size = 3, alpha = 1) +
  geom_line(size = 1, alpha = 1, aes(y=estimate, x = factor(week), color = treatment)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
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
  scale_x_discrete(name = "age (weeks)")


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

plot.morph.mm.mass <- ggplot() + 
  geom_point(size = 2.5, alpha = 0.7, data = morph.data.mm, aes(y=mass.g, x = days.forelimb, color = treatment)) +
  geom_ribbon(data = predicted.df.morph.mm,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group), show.legend = FALSE) +
  #scale_color_gradient(low = "gray85", high = "gray1") +
  geom_line(data = predicted.df.morph.mm, size = 1, aes(x = x, y = predicted, group = group, color = group)) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
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
  scale_y_continuous(name = "mass (g) at metamorphosis") +
  scale_x_continuous(name = "larval duration (days)", limits = c(42, 75), breaks = seq(40, 75, by = 5))

plot.morph.mm.smi.1 <- ggplot() + 
  geom_point(size = 3, alpha = 1, pch = 21, stroke = 1, data = morph.data.mm, aes(y=smi, x = mass.g, color = treatment, fill = days.forelimb)) +
 scale_fill_gradient(low = "gray85", high = "gray1") + 
 scale_color_manual(values = c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  geom_ribbon(data = predicted.df.bodycond.mass,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.bodycond.mass, size = 1, aes(x = x, y = predicted)) +
  theme_bw() +
  labs(color="larval duration (days)") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index at metamorphosis") +
  scale_x_continuous(name = "mass (g)")

plot.morph.mm.smi.2 <- ggplot() + 
  geom_jitter(width = 0.3, size = 3, stroke = 1, alpha = 1, pch = 21, data = morph.data.mm, aes(y=smi, x = days.forelimb, color = treatment, fill = mass.g)) +
  scale_fill_gradient(low = "gray85", high = "gray1") + 
  scale_color_manual(values = c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  geom_ribbon(data = predicted.df.bodycond.daysforelimb,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), show.legend = FALSE) +
  geom_line(data = predicted.df.bodycond.daysforelimb, size = 1, aes(x = x, y = predicted)) +
  theme_bw() +
  labs(fill="metamorphic mass (g)") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "scaled mass index at metamorphosis") +
  scale_x_continuous(name = "larval duration (days)", limits = c(42, 75), breaks = seq(40, 75, by = 5))

ggplot() + 
  geom_jitter(width = 0.3, size = 1, alpha = 0.6, pch = 21, stroke = 1, data = morph.data.mm, aes(y=smi, x = mean.days.forelimb, color = treatment, fill = treatment)) +
  geom_point(size = 3, alpha = 1, pch = 21, stroke = 1, data = morph.data.mm, color = "black", aes(y=mean.mm.smi, x = mean.days.forelimb, fill = treatment, group = juv.tank.id)) +
  scale_color_manual(values = c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values = c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  geom_ribbon(data = test,
              alpha = 0.4,
              mapping = aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group), show.legend = FALSE) +
  geom_line(data = test, size = 1, aes(x = x, y = predicted, color = group)) +
  theme_bw() +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c("left", "top"),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "mean tank scaled mass index at metamorphosis") +
  scale_x_continuous(name = "mean tank larval duration (days)")


# PLOT DATASETS: Effect of metamorphic mass on morphometrics at and after metamorphosis ---------------------

# x-y plot with observed and glmm-predicted values for mass
plot.morph.mass <- ggplot() + 
  # geom_ribbon(data = predicted.df.juv.morph,
  #             alpha = 0.75,
  #             mapping = aes(x = x, y = predicted, fill=factor(group, levels = c("large", "med", "small")), ymin = conf.low, max = conf.high), show.legend = FALSE) +
  geom_jitter(width = 0.8, size = 2, alpha = 1, data = morph.data.mm.juv, aes(y=mass.g, x = post.mm.weeks.num, color = mean.mm.mass.g, group = clutchtank)) + #individual level response colored by mean tank mass at metamorphosis
  # scale_fill_manual(values = c("gray1", "gray50", "gray85")) +
  geom_line(data = morph.data.mm.juv.summ, pch = 21, size = 0.8, aes(y=mean.mass.g, x = post.mm.weeks.num, color = mean.mm.mass.g, group = clutchtank)) +
  scale_color_gradient(low = "gray85", high = "gray1") +
  scale_fill_gradient(low = "gray85", high = "gray1") +
  
  new_scale_fill() +
  geom_point(data = morph.data.mm.juv.summ, pch = 21, alpha = 1, size = 4, stroke = 1, aes(y=mean.mass.g, x = post.mm.weeks.num, fill = treatment, color = mean.mm.mass.g,  group = clutchtank)) + #tank level response
  
  new_scale_color() +
  geom_line(data = predicted.df.morph.juv, color = "black", linetype = 2, size = 2, aes(x = x, y = predicted)) + #from model on individual mass but with log mass and poly3 week
  geom_line(data = predicted.df.morph.juv.summ, size = 2, linetype = 3, aes(x = x, y = predicted)) + #from model on tank-level means for mass but with log mass and poly3 week
  
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[3], natparks.pals("BryceCanyon")[1])) +

  scale_linetype_manual(values=c(1,2,3)) +
  theme_bw() +
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
  scale_y_continuous(name = "mass (g)") + 
  scale_x_continuous(name = "post-metamorphic age (weeks)", breaks = unique(morph.data.mm.juv$post.mm.weeks.num), labels = unique(morph.data.mm.juv$post.mm.weeks))





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


# PANEL PLOT DATASETS: Create panel plot with mass as a function of treatment and larval duration across all sampling points --------------
ggarrange(plot.morph1, plot.morph2, plot.morph3, plot.morph4,
          ncol = 2,
          nrow = 2,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))

ggarrange(plot.morph.tad.tl.ld, plot.morph.tad.tl.hd, 
          #plot.morph.tad.bl.ld, plot.morph.tad.bl.hd, 
          ncol = 2,
          nrow = 2,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))

ggarrange(plot.morph.tad.growth.tl.ld, plot.morph.tad.growth.tl.hd, 
          plot.morph.tad.growth.bl.ld, plot.morph.tad.growth.bl.hd, 
          ncol = 2,
          nrow = 2,
          common.legend = FALSE,
          legend = NULL,
          labels = c("a", "b", "c", "d"),
          font.label = list(size = 20, color = "black"))

ggarrange(plot.morph.mass, labels = c("a"), font.label = list(size = 20, color = "black"),
          ggarrange(plot.smi.1, plot.smi.2, labels = c("b", "c"),font.label = list(size = 20, color = "black", ncol = 2)),
          nrow = 2,
          common.legend = FALSE,
          legend = NULL
          )


# SUMMARY TABLES: Morphological metrics at and after metamorphosis ---------------------

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
  

