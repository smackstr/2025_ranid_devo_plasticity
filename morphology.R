# this script creates figures of the morphology results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on morphology metrics, 2. effects of larval duration on morphology metrics) across three life stages (1. larval period, 2. at metamorphosis, 3. up to 3 months post-metamorphosis)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)

setwd("~/Desktop/R Working Directory/Databases")


# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for metamorph morphometrics, juvenile morphometrics, and developmental timing data
morph.data.mm = read.csv("Database_Morphometrics - Metamorph Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
morph.data.juv = read.csv("Database_Morphometrics - Froglet_Toadlet Morphometrics.csv", header = TRUE, skip = 0, na.strings = "NA")
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
morph.data.mm = morph.data.mm[morph.data.mm$gs.code == "RS",]
morph.data.juv = morph.data.juv[morph.data.juv$gs.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS",]

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
morph.data.juv$post.mm.weeks = factor(morph.data.juv$post.mm.weeks, levels = c("1-2", "4-6", "5-7", "8-10", "11-12", "12-14"))

# create unique tank id for all datasets
morph.data.mm$unique.id = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$juv.tank.id, sep = "_")

# need to change juvenile tanks with "extra" in their name to be just the name of the tank
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J002extra"] = "J002"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J003extra"] = "J003"
morph.data.juv$juv.tank.id[morph.data.juv$juv.tank.id == "J005extra"] = "J005"
morph.data.juv$unique.id = paste(morph.data.juv$gs.code, morph.data.juv$clutch, morph.data.juv$juv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$juv.tank.id, sep = "_")

# create unique individual id for devo.data and morph.data.mm, since these are the two we have individual-level data for
morph.data.mm$unique.id.indiv = paste(morph.data.mm$gs.code, morph.data.mm$clutch, morph.data.mm$larv.tank.id, morph.data.mm$animal.id, sep = "_")
devo.data$unique.id.indiv = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, devo.data$animal.id, sep = "_")


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.juv -----------------------

# remove overflow individuals
morph.data.juv = morph.data.juv[morph.data.juv$treatment != "overflow",]

# create column to store the developmental data
morph.data.juv$mean.days.forelimb = NA

# fill developmental data column - average mean.days.forelimb for each juvenile tank
for(i in 1:length(unique(morph.data.juv$unique.id))){
  morph.data.juv$mean.days.forelimb[morph.data.juv$unique.id == unique(morph.data.juv$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(morph.data.juv$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to morph.data.mm -----------------------

# remove overflow individuals
morph.data.mm = morph.data.mm[morph.data.mm$treatment != "overflow",]

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
morph.data.mm.juv <- rbind(morph.data.mm[, c("date.measured",
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


# PLOT DATASETS: Effect of rearing density on morphometrics at and after metamorphosis -----------------------
# option 1 = x-y plot with summarized mean and +/- 1 se
ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=mass.g, x = post.mm.weeks, color = treatment)) + 
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
        axis.text.x=element_text(size=18, color = "black"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x=element_text(size=18, color = "black"), 
        axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "mass (g)") +
  scale_x_discrete(name = "weeks after metamorphosis")

ggplot(data = morph.data.mm.juv[!is.na(morph.data.mm.juv$juv.tank.id),], aes(y=mass.g, x = post.mm.weeks, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_boxplot(alpha = 0.75, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=18, color = "black"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x=element_text(size=18, color = "black"), 
        axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  expand_limits(y = 0) +
  scale_y_continuous(name = "mass (g)") +
  scale_x_discrete(name = "weeks after metamorphosis")



# PLOT DATASETS: Effect of developmental speed on morphometrics at and after metamorphosis ---------------------
