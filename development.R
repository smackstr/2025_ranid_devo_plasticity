# this script creates figures of the developmental results for the 2023 ranid developmental paper
# includes effects of rearing density on developmental timepoints (1. time to forelimb emergence, time to complete metamorphosis, laterality of forelimb emergence)

rm(list = ls())

library(dplyr)
library(ggpubr)
library(NatParksPalettes)
library(xfun)

setwd("~/Desktop/R Working Directory/Databases")

# COMPILE DATASETS: Prepare Datasets  -----------------------

# read in databases for developmental data
devo.data = read.csv("Database_Metamorphosis - Metamorphosis Log.csv", header = TRUE, skip = 0, na.strings = "NA")

# subset databases to only include Rana sylvatica
devo.data = devo.data[devo.data$gs.code == "RS",]

#change "control" to "low density" - called the low density treatment "CO for "control" during experiments because it has very different letters than "HD" for "high density", so was less likely to get confused. But technically it's a low density vs. high density comparison
devo.data$treatment[devo.data$treatment == "control"] = "low density"

#change column classes
devo.data$first.six = factor(devo.data$first.six, levels = c("yes", "no"))
devo.data$treatment = factor(devo.data$treatment)
devo.data$larv.tank.id = factor(devo.data$larv.tank.id)
devo.data$juv.tank.id = factor(devo.data$juv.tank.id)
devo.data$days.forelimb.tail = as.integer(devo.data$days.forelimb.tail)


# PLOT DATASETS: Effect of rearing density on days to developmental metrics at and after metamorphosis -----------------------

# plotted separately by clutch 
plot.devo.1 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = clutch, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_boxplot(alpha = 0.75, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=18, color = "black"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x=element_text(size=18, color = "black"), 
        axis.title.y = element_text(size=18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days)") +
  scale_x_discrete(name = "clutches separated")

# plotted with clutch lumped together under each density group 
plot.devo.2 <- ggplot(data = devo.data[devo.data$treatment != "overflow" & devo.data$first.six == "yes",], aes(y=days.forelimb, x = treatment, color = treatment)) + 
  geom_point(position=position_jitterdodge(), size = 2.5, alpha = 0.7) +
  geom_boxplot(alpha = 0.75, size = 0.75) +
  scale_color_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  scale_fill_manual(values=c(natparks.pals("BryceCanyon")[-2])) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=18, color = "white"), 
        axis.text.y=element_text(size=18, color = "black"), 
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_y_continuous(name = "larval duration (days)") +
  scale_x_discrete(name = "clutches combined")

ggarrange(plot.devo.1, plot.devo.2, 
          common.legend = TRUE, 
          legend = "bottom", 
          labels = c("a", "b"),
          font.label = list(size = 20, color = "black"))



# add the legend as the final plot within plot list so that it can be graphed within the grid
plotList.morph[[length(plotList.morph) + 1]] <- as_ggplot(get_legend(plotList.morph))

# create panel plot with all morphometrics data across all sampling points
ggarrange(plotlist = plotList.morph, legend = "none")