


#==============================================================================#
# Changes
# Includes external fixation duration and action trigger RT in dataset
#==============================================================================#

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)

# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

source("./Analysis/helper_functions.R")

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

# Remove unnecessary columns
data <- data %>%
  dplyr::select(ID, Acc, condition, block, orientation,
                foreperiod, RT, counterbalance, 
                extFixationDuration, action_trigger.rt)

# Coerce to factors
data <- data %>%
  mutate(across(c(ID, condition, block, orientation, foreperiod, counterbalance), as_factor))

data$condition <- data$condition %>%
  fct_relevel(c("external", "action"))

# Remove practice trials
data <- data %>%
  filter(condition != 'practice')

# Add trial number across the whole experiment and by block
data <- data %>%
  group_by(ID) %>%
  mutate(trial = seq(n())) %>%
  ungroup() %>%
  group_by(ID, block) %>%
  mutate(trial_bl = seq(n())) %>%
  ungroup()

# Create column for previous orientation and for comparison of current and previous orientations
data <- data %>%
  mutate(seqOri = ifelse(lag(orientation)==orientation, 'same', 'different'),
         prevOri = lag(orientation)) %>%
  mutate(seqOri = as.factor(seqOri),
         prevOri = as.factor(prevOri))

# Save data with error trials to assess accuracy
dataAcc <- data

# Keep only trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Coerce foreperiod back to numeric
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
dataAcc$numForeperiod <- as.numeric(as.character(dataAcc$foreperiod))

# Create variable for error rate
dataAcc$Error <- abs(dataAcc$Acc - 1)

# Create log10 of continuous independent variables
data$numLogFP <- log10(data$numForeperiod)
data$logFP <- as.factor(data$numLogFP)

dataAcc$numLogFP <- log10(dataAcc$numForeperiod)
dataAcc$logFP <- as.factor(dataAcc$numLogFP)

# Factor version of accuracy/error rate
dataAcc$acc_result <- as.factor(dataAcc$Acc)
dataAcc$error_result <- as.factor(dataAcc$Error)

# Remove extreme values
ntrials_before_extrem <- nrow(data)
data <- data %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)
ntrials_after_extrem <- nrow(data)

# Transform RT to reduce skew
data$logRT <- ifelse(!is.na(data$RT), log10(data$RT), NA) # log-transform
data$invRT <- ifelse(!is.na(data$RT), 1/data$RT, NA)

# Trimming
data2 <- data %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# No trimming
data <- data %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  ungroup()


# Average data
summaryData <- data %>%
  group_by(ID,foreperiod,condition,
           orientation,seqOri, prevOri,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            varRT = var(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

summaryData2 <- data2 %>%
  group_by(ID,foreperiod,condition,
           orientation, seqOri, prevOri,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            varRT = var(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

summaryDataAcc <- dataAcc %>%
  group_by(ID, foreperiod, condition) %>%
  summarise(meanAcc = mean(Acc),
            varAcc = var(RT),
            errorRate = mean(Error)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2)

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")

