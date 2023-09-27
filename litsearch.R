# this script creates a methodology figure for the 2023 ranid developmental paper
# compares our methodology with those of other published studies on ranids using rearing density as a treatment

rm(list = ls())

library(dplyr)
library(ggpubr)
library(ggstar)

setwd("~/Desktop/R Working Directory/Databases")


# COMPILE DATASETS: Prepare Datasets  -----------------------
data.lit = read.csv("Database_LitSearch - Sheet1.csv", header = TRUE, skip = 0)

# subset databases to only include Rana sylvatica
data.lit = data.lit[data.lit$family== "Ranidae",]

#change column classes
data.lit$pub_index = factor(data.lit$pub_index)


# PLOT DATASETS: Study design -----------------------

# plot with two groupings: study and species
ggplot(data = data.lit, aes(y=delta_density_perliter, x = delta_num_tad, shape = pub_index, fill = species, color = species, group=interaction(pub_index, species))) + 
  geom_line(position = position_jitterdodge(jitter.width = 0,
                                            jitter.height = 0)) +
  geom_point(position = position_jitterdodge(jitter.width = 0,
                                             jitter.height = 0), size = 5, alpha = 1) + 
  geom_star(data = data.lit[data.lit$pub_citation == "this study",], position = position_jitterdodge(jitter.width = 0, jitter.height = 0), size = 9, alpha = 1, color = "black", fill = "black", aes(starshape=pub_index)) +
  geom_line(data = data.lit[data.lit$pub_citation == "this study",], color = "black") +
  scale_fill_hue() +
  scale_color_hue() +
  scale_shape_manual(values = 0:20) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=12),
        axis.text.x=element_text(size=16, color = "black"), 
        axis.text.y=element_text(size=16, color = "black"), 
        axis.title.x=element_text(size=16, color = "black"), 
        axis.title.y = element_text(size=16),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  guides(shape = "none") + 
  scale_y_continuous(name = "delta density per liter", limits = c(0, 30)) +
  scale_x_continuous(name = "delta number of tadpoles")