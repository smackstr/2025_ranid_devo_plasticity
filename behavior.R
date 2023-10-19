# this script creates figures of the behavior results for the 2023 ranid developmental paper
# includes two main branches (1. effects of rearing density on behavior metrics, 2. effects of larval duration on behavior metrics) across larval period

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)

setwd("~/Desktop/R Working Directory")
# read in Rdata output that has compiled all data from EthoVision outputs
load("EthoVision_CompileData.RData")

# COMPILE DATASETS: Prepare Datasets  -----------------------
setwd("~/Desktop/R Working Directory/Databases")

#read in metamorphosis timing log database
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA") 

# subset databases to only include Rana sylvatica
ev.data.summary.30 = ev.data.summary.30[ev.data.summary.30$genus.species.code == "RS",]
devo.data = devo.data[devo.data$gs.code == "RS",]

# change "control" to "low density"
ev.data.summary.30$treatment = as.character(ev.data.summary.30$treatment)
ev.data.summary.30$treatment[ev.data.summary.30$treatment == "control"] = "low density"
ev.data.summary.30$treatment = factor(ev.data.summary.30$treatment)

pres.abs.summ$treatment = as.character(pres.abs.summ$treatment)
pres.abs.summ$treatment[pres.abs.summ$treatment == "control"] = "low density"
pres.abs.summ$treatment = factor(pres.abs.summ$treatment)

devo.data$treatment[devo.data$treatment == "control"] = "low density"

# change column classes
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)

# subset devo dataset to only include first six individuals
devo.data = devo.data[devo.data$first.six == "yes",]

#temporarily removed NAs (need to figure out why they exist)
ev.data.summary.30 = ev.data.summary.30[is.na(ev.data.summary.30$stim.time.cat.30) == FALSE,]

# create unique tank id for ev.data.summary.30, pres.abs.summ, and devo.data
ev.data.summary.30$individual.id = as.character(ev.data.summary.30$individual.id)
ev.data.summary.30$larv.tank.id = NA
for(i in 1:nrow(ev.data.summary.30)){
  ev.data.summary.30$larv.tank.id[i] = strsplit(ev.data.summary.30$individual.id[i], split = "_")[[1]][3]
}
ev.data.summary.30$unique.id = paste(ev.data.summary.30$genus.species.code, ev.data.summary.30$clutch, ev.data.summary.30$larv.tank.id, sep = "_")

devo.data$unique.id = paste(devo.data$gs.code, devo.data$clutch, devo.data$larv.tank.id, sep = "_")

pres.abs.summ$individual.id = as.character(pres.abs.summ$individual.id)
pres.abs.summ$genus.species.code = NA
pres.abs.summ$clutch= NA
pres.abs.summ$larv.tank.id = NA
for(i in 1:nrow(pres.abs.summ)){
  pres.abs.summ$genus.species.code[i] = strsplit(pres.abs.summ$individual.id[i], split = "_")[[1]][1]
  pres.abs.summ$clutch[i] = strsplit(pres.abs.summ$individual.id[i], split = "_")[[1]][2]
  pres.abs.summ$larv.tank.id[i] = strsplit(pres.abs.summ$individual.id[i], split = "_")[[1]][3]
}
pres.abs.summ$unique.id = paste(pres.abs.summ$genus.species.code, pres.abs.summ$clutch, pres.abs.summ$larv.tank.id, sep = "_")


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to ev.data.summary.30 -----------------------

# create column to store the developmental data, mean.days.forelimb represents an average, since we don't have individual level data for this
ev.data.summary.30$mean.days.forelimb = NA

# fill developmental data column for survi.data.tad.exp
for(i in 1:length(unique(ev.data.summary.30$unique.id))){
  ev.data.summary.30$mean.days.forelimb[ev.data.summary.30$unique.id == unique(ev.data.summary.30$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(ev.data.summary.30$unique.id)[i]], na.rm = TRUE)
}


# COMPILE DATASETS: Add Developmental Data (mean days to forelimb) to pres.abs.summ -----------------------

# create column to store the developmental data, mean.days.forelimb represents an average, since we don't have individual level data for this
pres.abs.summ$mean.days.forelimb = NA

# fill developmental data column for survi.data.tad.exp
for(i in 1:length(unique(pres.abs.summ$unique.id))){
  pres.abs.summ$mean.days.forelimb[pres.abs.summ$unique.id == unique(pres.abs.summ$unique.id)[i]] <- mean(devo.data$days.forelimb[devo.data$unique.id ==  unique(pres.abs.summ$unique.id)[i]], na.rm = TRUE)
}



# PLOT DATASETS: Effect of rearing density on % frames in motion across three sampling points during larval period -----------------------

# option 1 = plotting only open field start
plot.behav.1 <- ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=100*(movement/frames), x = week, color = treatment)) + 
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
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,100)) +
  scale_x_discrete(name = "age (weeks)")


# option 2 = plotting across all stimulus time categories
ggplot(data = ev.data.summary.30, aes(y=100*(movement/frames), x = stim.time.cat.30, color = treatment)) + 
  facet_grid(rows = vars(week)) +
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
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,100)) +
  scale_x_discrete(name = "", labels = c("open field", "prelight", "light on", "postlight", "open field", "pretap", "posttap", "open field"))


# PLOT DATASETS: Effect of rearing density on % arena visited across three sampling points during larval period -----------------------

plot.behav.2 <- ggplot(data = pres.abs.summ, aes(y=100*prop.arena.visited, x = week, color = treatment)) + 
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
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena visited", limits = c(0,100)) +
  scale_x_discrete(name = "age (weeks)")

# PLOT DATASETS: Combine plots -----------------------
ggarrange(plot.behav.1, plot.behav.2, 
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))


# PLOT DATASETS: Effect of developmental speed on % frames in motion across three sampling points during larval period -----------------------

# option 1 = plotting only open field start
plot.behav.3 <- ggplot(data = ev.data.summary.30[ev.data.summary.30$stim.time.cat.30 == "open-field start",], aes(y=100*(movement/frames), x = mean.days.forelimb, color = week)) + 
  #facet_grid(rows = vars(week)) +
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% frames in motion", limits = c(0,100)) +
  scale_x_continuous(name = "larval duration (days)")

plot.behav.4 <- ggplot(data = pres.abs.summ, aes(y=100*prop.arena.visited, x = mean.days.forelimb, color = week)) + 
  #facet_grid(rows = vars(week)) +
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) + #fit linear model with confidence interval
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-1], natparks.pals("BryceCanyon")[-3])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x=element_text(size=16, color = "black", angle=0, hjust = 0.5, vjust = 0), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "% arena visited", limits = c(0,100)) +
  scale_x_continuous(name = "larval duration (days)")


# PLOT DATASETS: Combine plots -----------------------
ggarrange(plot.behav.3, plot.behav.4, 
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))