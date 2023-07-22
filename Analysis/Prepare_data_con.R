


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
                extFixationDuration, action_trigger.rt) %>%
  mutate(foreperiod = foreperiod * 1000,
         RT = RT *1000,
         extFixationDuration = extFixationDuration * 1000,
         action_trigger.rt = action_trigger.rt * 1000)
  

# Coerce to factors
data$ID <- as.factor(data$ID)
data$condition <- data$condition %>%
  as.factor() %>%
  forcats::fct_relevel(c("external", "action"))
data$block <- as.factor(data$block)
data$orientation <- as.factor(data$orientation)
data$foreperiod <- as.factor(data$foreperiod)
data$counterbalance <- as.factor(data$counterbalance)


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
dataAll <- data

# Keep only trials with correct responses to analyze RT
data <- data %>%
  filter(!is.na(RT), Acc == 1)

# Coerce foreperiod aback to numeric
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
dataAll$numForeperiod <- as.numeric(as.character(dataAll$foreperiod))

# Create log10 of continuous independent variables
data$numLogFP <- log10(data$numForeperiod)
data$logFP <- as.factor(data$numLogFP)

dataAll$numLogFP <- log10(dataAll$numForeperiod)
dataAll$logFP <- as.factor(dataAll$numLogFP)

# Factor version of accuracy
dataAll$acc_result <- as.factor(dataAll$Acc)

# Remove extreme values
data <- data %>%
  filter(RT < 1000) %>%
  filter(RT > 150)

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
  #filter(abs(logRTzscore) < 3) %>%
  ungroup()


# Average data
summaryData <- data %>%
  group_by(ID,foreperiod,condition,
           #orientation,seqOri, prevOri,
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
           #orientation, seqOri, prevOri,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            varRT = var(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

summaryDataAll <- dataAll %>%
  group_by(ID, foreperiod, condition) %>%
  summarise(meanAcc = mean(Acc),
            varAcc = var(RT)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2)

#write_csv(data, "./Analysis/data.csv")
#write_csv(data2, "./Analysis/data2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")

