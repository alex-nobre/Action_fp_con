
# Load necessary packages
library(tidyverse)
library(broom)
library(magrittr)
library(lattice)
library(afex)
library(emmeans)
library(lme4)
library(car)
library(data.table)
library(codingMatrices)
library(performance)
library(modelr)
library(BayesFactor)
library(bayestestR)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

# Load data
source('./Analysis/Prepare_data_con.R')


# Function to sample trials
trsample <- function(dataset, ntrials) {
  sampDS <- dataset %>%
    group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
    slice_sample(n = ntrials, replace = TRUE) %>% # Sample trials
    ungroup()
  
  summarySampDS <- sampDS %>%
    group_by(ID, foreperiod, condition) %>%
    summarise(meanRT = mean(RT)) %>%
    ungroup()
    
  sampDSAnova <- aov_ez(data = summarySampDS,
                        id = "ID",
                        dv = "meanRT", within = c("foreperiod", "condition"),
                        anova_table = list(es = "pes"))
  DSes <- sampDSAnova$anova_table[2,5]
}

# Create dataset with 50 sampled trials
sampdata <- data2 %>%
  group_by(ID, foreperiod, condition) %>% # Group by ID, FP, condition
  slice_sample(n = 12) %>% # Sample trials
  ungroup()

summarySampData <- sampdata %>%
  group_by(ID,foreperiod,condition) %>%
  summarise(meanRT = mean(RT)) %>%
  ungroup()

# analyze
ggplot(data = summarySampData, aes(x = foreperiod,
                            y = meanRT,
                            color = condition)) +
  stat_summary(fun = "mean", geom = "point",) +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c("orange", "blue"))

sampAnova = aov_ez(data = summarySampData,
                   id = "ID",
                   dv = "meanRT", within = c("foreperiod", "condition"),
                   anova_table = list(es = "pes"))


nSim <- 1000
esList <- vector(mode = "numeric", length = nSim)

options(dplyr.summarise.inform = FALSE)
for(iSim in 1:nSim) {
  esList[iSim] <- trsample(data2, 12)
}
options(options_defaults)

plot(1:nSim, esList)

mean(esList)


ntrList <- 8:14
nSim <- 1000
esByTR <- vector(mode = "list", length = length(ntrList))

options(dplyr.summarise.inform = FALSE)
for(thisTr in 1:length(esByTR)) {
  ntrials <- ntrList[thisTr]
  esList <- vector(mode = "numeric", length = nSim)
  for(iSim in 1:nSim) {
    esList[iSim] <- trsample(data2, ntrials)
  }
  esByTR[[thisTr]] <- esList
}

options(options_defaults)

par(mfrow = c(4,2))

jpeg("./Analysis/Plots/trials_to_effect.jpeg", width = 1200, height = 1000)
for(thisTr in 1:length(ntrList)) {
  esToPlot <- esByTR[[thisTr]]
  plot(1:nSim, esToPlot)
}
dev.off()
par(graphical_defaults)
