# README -------------------------
# this script (1) compares our methodology with those of other published studies on ranids using rearing density as a treatment, (2) creates a methodology figure for the 2025 ranid developmental paper

# LOAD PACKAGES and SET WD ------------------------------------
rm(list = ls())

library(dplyr)
library(ggpubr)
library(ggstar)
library(lme4)
library(car)

setwd("~/Desktop/R Working Directory/Databases")


# COMPILE DATASETS: Prepare Datasets  -----------------------
data.lit = read.csv("Database_LitSearch - Sheet1.csv", header = TRUE, skip = 0)

# subset databases to only include Rana sylvatica
data.lit = data.lit[data.lit$family== "Ranidae",]

#change column classes
data.lit$pub_index = factor(data.lit$pub_index)

#create column to say whether sylvatica or not
data.lit$sylvatica[data.lit$species == "sylvatica"] = "sylvatica"
data.lit$sylvatica[data.lit$species != "sylvatica"] = "other ranid"
data.lit$sylvatica = factor(data.lit$sylvatica)

#specify factor so that colors assigned correctly in plots
data.lit$species = factor(data.lit$species, levels = c("sylvatica", "temporaria", "pipiens", "lessonae", "lessonae x ridibundus", "ridibundus", "sphenocephala", "temporalis", "areolata"))

#create shortened version of pub_citation
data.lit$pub_citation = str_extract(data.lit$pub_citation, "[^([[:digit:]]).]+")
sub("([[:digit:]])._.*", "", data.lit$pub_citation)


#ANALYZE DATASETS: Does change in response scale with severity -----------------------
lmer.full <- lmer(delta_mass ~ scale(delta_density_perliter)*species + (1|pub_citation), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.nointxn <- lmer(delta_mass ~ scale(delta_density_perliter) + species + (1|pub_citation), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.nospecies <- lmer(delta_mass ~ scale(delta_density_perliter) + (1|pub_citation), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.noseverity <- lmer(delta_mass ~ species + (1|pub_citation), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.null <- lmer(delta_mass ~ (1|pub_citation), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

anova(lmer.full, lmer.nointxn, lmer.nospecies, lmer.noseverity, lmer.null)

summary(lmer.nointxn)
Anova(lmer.nointxn, type = "II")

lmer.nointxn <- lmer(delta_larval_dur ~ scale(delta_density_perliter) + species_reported + (1|pub_index), data = data.lit[data.lit$family == "Ranidae" & data.lit$delta_density_perliter > 0,], na.action = na.omit)
summary(lmer.nointxn)
Anova(lmer.nointxn, type = "II")



#ANALYZE DATASETS: Does change in response scale with severity ONLY FOR WOOD FROGS -----------------------
lmer.full <- lmer(delta_mass ~ delta_density_perliter + (1|pub_index), data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.null <- lmer(delta_mass ~ (1|pub_index), data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

anova(lmer.full, lmer.null)

summary(lmer.full)
confint(lmer.full)


lmer.full <- lmer(delta_larval_dur ~ delta_density_perliter + (1|pub_citation), data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

lmer.null <- lmer(delta_larval_dur ~ (1|pub_citation), data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_density_perliter > 0,], na.action = na.omit)

anova(lmer.full, lmer.null)

summary(lmer.full)
confint(lmer.full)



# PLOT DATASETS: Study design -----------------------

# plot with two groupings: study and species
plot2 <- ggplot(data = data.lit[data.lit$pub_citation != "this study",], aes(y=delta_density_perliter, x = delta_num_tad, group=interaction(pub_index, sylvatica))) + 
  geom_line(position = position_dodge(5), aes(color = sylvatica)) +
  geom_point(size = 2, alpha = 0.4, shape = 21, color = "black", position = position_dodge(5), aes(fill = sylvatica)) + 
  geom_star(alpha = 0.4, data = data.lit[data.lit$pub_citation == "this study",], size = 5, color = "black", fill = "cornflowerblue", aes(starshape=pub_index)) +
  geom_line(data = data.lit[data.lit$pub_citation == "this study",], color = "cornflowerblue") +
  scale_color_manual(values = c("gray45", "cornflowerblue")) +
  scale_fill_manual(values = c("gray45", "cornflowerblue")) +
  #scale_color_hue() +
  #scale_shape_manual(values = 0:20) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = "number of tadpoles", limits = c(0,300))


# without lines ---------

ggplot(data = data.lit %>% filter(pub_citation != "this study",
                                  delta_num_tad > 0), 
       aes(y=delta_density_perliter, x = delta_num_tad, group=interaction(pub_index, sylvatica))) + 

  geom_point(size = 2, alpha = 0.4, shape = 21, color = "black", position = position_dodge(5), aes(fill = sylvatica)) + 
 
  #our study
   geom_star(alpha = 1, data = data.lit %>% filter(pub_citation == "this study",
                                                     delta_num_tad > 0), 
             size = 5, color = "black", fill = "cornflowerblue", aes(starshape=pub_index)) +

  geom_vline(xintercept = 50, linetype="dotted", 
             color = "black", size=0.8) +
  
  geom_hline(yintercept = 10, linetype="dotted", 
             color = "black", size=0.8) +
  
  scale_color_manual(values = c("gray45", "cornflowerblue")) +
  scale_fill_manual(values = c("gray45", "cornflowerblue")) +
  #scale_color_hue() +
  #scale_shape_manual(values = 0:20) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "number of tadpoles", sep=" "), limits = c(0,250))



#only in wood frogs ---------
ggplot(data = data.lit[data.lit$pub_citation != "this study" & data.lit$species == "sylvatica",], aes(y=delta_density_perliter, x = delta_num_tad, group=interaction(pub_index, sylvatica))) + 
  geom_line(position = position_dodge(5), aes(color = sylvatica)) +
  geom_point(size = 2, alpha = 0.4, shape = 21, color = "black", position = position_dodge(5), aes(fill = sylvatica)) + 
  geom_star(alpha = 0.4, data = data.lit[data.lit$pub_citation == "this study",], size = 5, color = "black", fill = "cornflowerblue", aes(starshape=pub_index)) +
  geom_line(data = data.lit[data.lit$pub_citation == "this study",], color = "cornflowerblue") +
  scale_color_manual(values = c("gray45", "cornflowerblue")) +
  scale_fill_manual(values = c("gray45", "cornflowerblue")) +
  #scale_color_hue() +
  #scale_shape_manual(values = 0:20) +
  theme_bw() +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.title = element_blank(),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = "number of tadpoles", limits = c(0,300))


# PLOT DATASETS: Study results within context of study design: ranids vs. other ranids (not separated by species) -----------------------

# results for change in mass
plot3 <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # other ranid lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$species != "sylvatica" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, size = num_tad, group=interaction(pub_index, species)),
             color = "gray45") + 
  
  #sylvatica lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, size = num_tad, group=interaction(pub_index, species)),
             color = "cornflowerblue") + 
  
  #our study point + star
  geom_star(alpha = 1, data = data.lit[data.lit$pub_citation == "this study"& data.lit$delta_mass != 0,], size = 10.5, color = "cornflowerblue", fill = "cornflowerblue", aes(starshape=pub_index, y=delta_density_perliter, x = delta_mass), show.legend = FALSE) +

  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_mass != 0,], alpha = 1, color = "black", fill = "cornflowerblue", shape = 21, aes(y=delta_density_perliter, x = delta_mass, size = num_tad), show.legend = FALSE) + 

  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 320), 
                        breaks = seq(5,320,45)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none", size=guide_legend(title="number of tadpoles")) + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "mass", sep=" "), limits = c(-1.5,1.0), breaks = seq(-1.5,1.0,0.25))


#results for change in larval duration
plot4 <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # other ranid lines
  #geom_line(data = data.lit[data.lit$species != "sylvatica" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, group=interaction(pub_index, species)),
  #color = "gray45") +
  geom_point(alpha = 0.6, data = data.lit[data.lit$species != "sylvatica" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad, group=interaction(pub_index, species)),
             color = "gray45") + 
  
  #sylvatica lines
  #geom_line(data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, group=interaction(pub_index, species)),
  #color = "cornflowerblue") +
  geom_point(alpha = 0.6, data = data.lit[data.lit$species == "sylvatica" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad, group=interaction(pub_index, species)),
             color = "cornflowerblue") + 
  
  #our study point + star
  geom_star(alpha = 1, data = data.lit[data.lit$pub_citation == "this study"& data.lit$delta_mass != 0,], size = 10.5, color = "cornflowerblue", fill = "cornflowerblue", aes(starshape=pub_index, y=delta_density_perliter, x = delta_larval_dur), show.legend = FALSE) +
  
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0,], alpha = 1, color = "black", fill = "cornflowerblue", shape = 21, aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 320), 
                        breaks = seq(5,320,45)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none", size=guide_legend(title="number of tadpoles")) + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "larval duration", sep=" "), limits = c(-15,45), breaks = seq(-15,45,5))

ggarrange(plot3,plot4,
          common.legend = TRUE,
          legend = "bottom",
          labels = c("a", "b"))


# PLOT DATASETS: Study results within context of study design: ranids vs. other ranids (separated by species) -----------------------

# results for change in mass
plot3 <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # all lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, size = num_tad, color = species, group=interaction(pub_index, species))) +

  geom_text(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = delta_mass, label = pub_index)) +
  
  #our study line + star
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_mass != 0,], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = delta_mass, size = num_tad), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "mass", sep=" "), limits = c(-1.5,0.5), breaks = seq(-1.5,0.5,0.25))


#results for change in larval duration
plot4 <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # other ranid lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad, color = species, group=interaction(pub_index, species))) + 
  
  geom_text(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = delta_larval_dur, label = pub_index)) + 
  
  #our study line + star
  #geom_line(data = data.lit[data.lit$pub_citation == "this study",], aes(y=delta_density_perliter, x = num_tad),
  #color = "cornflowerblue") +
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0,], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "larval duration", sep=" "), limits = c(-10,30), breaks = seq(-10,30,5))


# just with wood frogs
rs.plot.days <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # add colored points by pub citation without showing scale
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0 & data.lit$sylvatica == "sylvatica",], 
             pch=21,
             aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad, color = pub_index, fill = pub_index), show.legend = FALSE) + 
  
  # add outlines points by pub citation 
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0 & data.lit$sylvatica == "sylvatica",], 
             color = "black", pch=21, stroke = 1, 
             aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad), show.legend = TRUE) + 
  
  
  #our study line
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0 & data.lit$sylvatica == "sylvatica",], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad), show.legend = FALSE) + 
  
  #scale_fill_hue(drop=TRUE) +
  #scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25),
                        name = "number of tadpoles\nin high density") +
  theme_bw() +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        axis.text.x=element_text(size=18, color = "black"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x=element_text(size=18, color = "black"), 
        axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none", fill = FALSE, color = FALSE) + #remove legend for shape and fill
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "larval duration", sep=" "), limits = c(-10,30), breaks = seq(-10,30,5))



rs.plot.days.LEGEND <- ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # add colored points by pub citation without showing scale
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0 & data.lit$sylvatica == "sylvatica",], 
             pch=21,
             aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad, fill = pub_citation), show.legend = TRUE) + 
  
  #our study line
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0 & data.lit$sylvatica == "sylvatica",], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = delta_larval_dur, size = num_tad), show.legend = FALSE) + 
  
  #scale_fill_hue(drop=TRUE) +
  #scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25),
                        name = "number of tadpoles\nin high density",
                        labels = function(x) str_wrap(x, width = 2)) +
  theme_bw() +
  theme(legend.title = element_text(size=16),
        legend.key.height = unit(2, "cm"),
        legend.margin = margin(),
        legend.text = element_text(size=8),
        axis.text.x=element_text(size=18, color = "black"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x=element_text(size=18, color = "black"), 
        axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none", size = FALSE) + #remove legend for shape and size
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("\u0394", "larval duration", sep=" "), limits = c(-10,30), breaks = seq(-10,30,5))


legend <- cowplot::get_legend(rs.plot.days)
grid::grid.draw(legend)

# PANEL PLOT DATASETS: create panel plot with relative and absolute tadpole perspectives
ggarrange(plot2, legend = NULL, labels = c("a"), font.label = list(size = 20, color = "black"),
          ggarrange(plot3, plot4, common.legend = TRUE, legend = "bottom", labels = c("b", "c"), font.label = list(size = 20, color = "black")),
          nrow = 2
          )

ggarrange(plot3, plot4, common.legend = TRUE, legend = "bottom", labels = c("b", "c"), font.label = list(size = 20, color = "black"),
          nrow = 1)

# PLOT DATASETS: LOG ODDS within context of study design: ranids vs. other ranids (separated by species) -----------------------

# results for change in mass
ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # all lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_mass != 0,], aes(y=delta_density_perliter, x = log_odds_mass, size = num_tad, color = species, group=interaction(pub_index, species))) +
  
  geom_text(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_mass != 0,],aes(y=delta_density_perliter, x = log_odds_mass, label = pub_index)) +
  
  #our study line + star
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_mass != 0,], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = log_odds_mass, size = num_tad), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("log odds mass", sep=" "))


#results for change in larval duration
ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # other ranid lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = log_odds_duration, size = num_tad, color = species, group=interaction(pub_index, species))) + 
  
  geom_text(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], aes(y=delta_density_perliter, x = log_odds_duration, label = pub_index)) + 
  
  #our study line + star
  #geom_line(data = data.lit[data.lit$pub_citation == "this study",], aes(y=delta_density_perliter, x = num_tad),
  #color = "cornflowerblue") +
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0,], alpha = 1, color = "black", fill = "black", shape = 21, aes(y=delta_density_perliter, x = log_odds_duration, size = num_tad), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = paste("\u0394", "density per liter", sep=" "), limits = c(0, 30)) +
  scale_x_continuous(name = paste("log odds larval duration", sep=" "))


#results for change in larval duration AND MASS, grouped by study
ggplot() + 
  #vertical line at 0 (no change in response variable)
  geom_vline(xintercept = 0, linetype="dashed", 
             color = "black", size=0.5) +
  
  # other ranid lines
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], shape = 2,
             aes(y=species, x = log_odds_duration, color = species, group=species)) + 
  
  geom_point(alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], shape = 1,
             aes(y=species, x = log_odds_mass, color = species, group=species)) + 
  
  # mean for other ranid lines
  stat_summary(fun.y=mean, alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], shape = 24, color = "black", fill = "black",
             aes(y=species, x = log_odds_duration, group=species)) + 
  
  stat_summary(fun.y=mean, alpha = 0.6, data = data.lit[data.lit$pub_citation != "this study" & data.lit$delta_larval_dur != 0,], shape = 21, color = "black", fill = "black",
             aes(y=species, x = log_odds_mass, group=species)) + 

  
  #our study line + star
  #geom_line(data = data.lit[data.lit$pub_citation == "this study",], aes(y=delta_density_perliter, x = num_tad),
  #color = "cornflowerblue") +
  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0,], alpha = 1, color = "black", fill = "black", shape = 2,
             aes(y=species, x = log_odds_duration, group = species, show.legend = FALSE)) + 

  geom_point(data = data.lit[data.lit$pub_citation == "this study" & data.lit$delta_larval_dur != 0,], alpha = 1, color = "black", fill = "black", shape = 1,
             aes(y=species, x = log_odds_mass, group = species), show.legend = FALSE) + 
  
  scale_fill_hue() +
  scale_color_hue() +
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 250), 
                        breaks = seq(5,255,25)) +
  theme_bw() +
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        axis.text.x=element_text(size=12, color = "black"), 
        axis.text.y=element_text(size=12, color = "black"), 
        axis.title.x=element_text(size=12, color = "black"), 
        axis.title.y = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none")
