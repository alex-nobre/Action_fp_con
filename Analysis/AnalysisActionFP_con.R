

#================================================================================#
# Changes
# Manually set contrasts for anova
#================================================================================#

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

# Load data
source('./Analysis/Prepare_data_con.R')

#==========================================================================================#
#======================================= 0. Data quality ===================================
#==========================================================================================#
# 1.1.1. RTs across blocks (conditions aggregated)
ggplot(data=summaryData2,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line', linewidth = 1, aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  labs(title='RT by block')

# 1.1.2. RTs across blocks (separated by condition)
ggplot(data=summaryData2,
       aes(x=block,
           y=meanRT,
           color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',linewidth=1,aes(group=condition))+
  stat_summary(fun.data='mean_sdl',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('orange','blue')) +
  labs(title='RT by block and condition')

# 1.1.3. RTs across blocks by counterbalancing order
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT,
           color=counterbalance))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean', geom='line', linewidth = 1, aes(group=counterbalance))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('blue','orange')) +
  labs(title='RT by block split by counterbalancing order')

# Distribution of data
dataHists <- ggplot(data=summaryData,
                    aes(x=meanRT,
                        color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

dataHists

# Individual histograms
indHistograms <- ggplot(data=data,
                        aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
indHistograms

qqmath(~RT|ID, data=data)
qqmath(~invRT|ID, data=data)

# 1.4.1. Plot RT by foreperiod only
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line' ) +
  stat_summary(fun.data = 'mean_cl_boot', width = 0.2, geom = 'errorbar') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~ID)

# Check for influence of external fixation duration
ggplot(data=filter(data,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod+ID)

# Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)

#================================ 0.2. Stopping rule =========================================
#======= 0.2.1. Plots as a function of sample size ======
plotsList <- list()

for(part in 2:length(unique(data2$ID))) {
  parts <- unique(data2$ID)[1:part]
  plotData <- summaryData2 %>%
    filter(ID %in% parts)
  
  if(part == length(unique(data2$ID))){
    thisPlot <- ggplot(data = plotData,
                       aes(x = foreperiod,
                           y = meanRT,
                           color = condition)) +
      stat_summary(fun = "mean", geom = "point") +
      stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
      stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
      labs(title = paste(c(part, "parts"), collapse = " "),
           x = "Foreperiod",
           y = "Mean RT") +
      scale_color_manual(values = c("orange","blue"))
  } else {
    thisPlot <- ggplot(data = plotData,
                       aes(x = foreperiod,
                           y = meanRT,
                           color = condition)) +
      stat_summary(fun = "mean", geom = "point") +
      stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
      stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
      labs(title = paste(c(part, "parts"), collapse = " "),
           x = "Foreperiod",
           y = "Mean RT") +
      scale_color_manual(values = c("orange","blue")) +
      theme(legend.position = "none") 
  }
  plotsList[[part-1]] <- thisPlot
}

stopRulePlot <- cowplot::plot_grid(plotlist = plotsList)
ggsave("./Analysis/Plots/seq_analysis.png",
       stopRulePlot,
       width = 15.0,
       height = 9.0)

#========= 0.2.3. Individual linear models comparison =========================

# Variables used as predictors: numForeperiod
# Dependent variable: RT
# Variables nested by condition and ID

buildmodel <- function(data, RT) {
  lm(RT ~ numForeperiod,
     data = data)
}

# Nest data to fit models
nested_data <- data2 %>%
  select(ID, condition, numForeperiod, RT) %>%
  group_by(ID, condition) %>%
  nest()

# Fit models and save parameters as columns
fitted_data <- nested_data %>%
  mutate(fit = map(data, buildmodel),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, condition, term, estimate) %>%
  pivot_wider(names_from = term,
              values_from = estimate)


# Foreperiod
fp_bfs <- ttestBF(x = fitted_data$numForeperiod[fitted_data$condition=='external'],
                  y = fitted_data$numForeperiod[fitted_data$condition=='action'],
                  paired=TRUE)

#========= 0.2.3. Mixed models BF comparison ============
library(buildmer)

with_condition <- buildmer(RT ~ numForeperiod * condition + (1 + numForeperiod * condition | ID),
                           data = data2,
                           buildmerControl = list(direction = "backward",
                                                  crit = "LRT",
                                                  family = gaussian(link = "identity"),
                                                  calc.anova = TRUE))

isSingular(with_condition)

formula(with_condition)

with_condition <- mixed(formula = RT ~ 1 + condition + foreperiod + condition:foreperiod + 
                          (1 + condition + foreperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = RT ~ 1 + foreperiod + 
                        (1 + foreperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

#========= 0.2.4. Sequential bayes factor ============
external_fits <- fitted_data[fitted_data$condition=='external',]
action_fits <- fitted_data[fitted_data$condition=='action',]

srange <- 10:nrow(external_fits)

fp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$numForeperiod[1:range],
                    y = action_fits$numForeperiod[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

jpeg("./Analysis/Plots/BF_seq.jpeg", width = 600, height = 500)
plot(srange, fp_bfs,
     xlab = "n of participants",
     ylab = "Bayes Factor",
     main = "Bayes Factors for FP effect compared between conditions")
lines(srange, fp_bfs)
dev.off()

#==========================================================================================#
#==================================== 1. Descriptives ======================================
#==========================================================================================#

meandata <- summaryData2 %>%
  group_by(condition, foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandata  

#==========================================================================================#
#==================================== 2. Basic models ======================================
#==========================================================================================#

# Set constrasts for variables used in ANOVAs
contrasts(summaryData$foreperiod) <- c(-1/2, 1/2)
contrasts(summaryData$condition) <- c(-1/2, 1/2)

contrasts(summaryData2$foreperiod) <- c(-1/2, 1/2)
contrasts(summaryData2$condition) <- c(-1/2, 1/2)

#==================== 2.1. FP x RT by condition ======================
# Lines by condition
lines_by_condition <- ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange","blue")) +
  facet_wrap(~ID)
ggplot2::ggsave("./Analysis/Plots/plot_by_sub.png",
                lines_by_condition)

lines_by_condition <- ggplot(data = summaryData2,
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange","blue"))
ggplot2::ggsave("./Analysis/Plots/plot_by_condition.png",
                lines_by_condition,
                width = 7.7,
                height = 5.8)


# Plot error bars based on differences from subject average
condDiffs <- function(column) {
  # Compute sub RT mean
  mRT <- mean(RT)
}

# diffsData <- summaryData2 %>%
#   group_by(ID) %>%
#   summarise(globalRT = mean(meanRT))

diffsData <- summaryData2 %>%
  group_by(ID) %>%
  plyr::ddply("ID", transform, grMeanRT = mean(meanRT)) %>%
  group_by(ID, condition, foreperiod) %>%
  summarise(meanRT = mean(meanRT), grMeanRT = mean(grMeanRT)) %>%
  ungroup() %>%
  mutate(condsRT = meanRT - grMeanRT)

ggplot(data = diffsData,
       aes(x = foreperiod,
           y = condsRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", width = 0.2, geom = "errorbar") +
  #stat_summary(fun.data = "mean_sd", width = 0.2, geom = "errorbar") +
  labs(title = "RT differences from mean (se as errorbars)",
       x = "Foreperiod",
       y = "RT") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange","blue"))
  

 # Anova
fpAnova <- aov_ez(id = "ID",
       dv = "meanRT",
       data = summaryData2,
       within = c("foreperiod", "condition"),
       anova_table = list(es = "pes"))

### Check assumptions

# Sphericity
check_sphericity(fpAnova)

fpAnova_plot <- afex_plot(fpAnova, x = 'foreperiod', trace = 'condition', error = 'within')

# Normality of residuals
is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")
plot(is_norm, type = "qq", detrend = TRUE)

testnormality = function(dfr) return(shapiro.test(dfr$invRT)$p.value)
p = as.vector(by(data, data$ID, testnormality))
names(p) = levels(data$ID)
names(p[p < 0.05])


# Try transformations
invfpAnova <- aov_ez(id = "ID",
                  dv = "meanInvRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(invfpAnova)

# Normality of residuals
is_norm <- check_normality(invfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# using 1/RT does not solve the problem

logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))

### Check assumptions

# Sphericity
check_sphericity(logfpAnova)

# Normality of residuals
is_norm <- check_normality(logfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# The log-transform does not solve the problem either

fpregression <- lm(meanRT ~ condition * foreperiod, data = summaryData)
summary(fpregression)
anova(fpregression)

logfpregression <- lm(meanRT ~ condition * logFP, data = summaryData)
anova(logfpregression)

# fpEmmeans <- emmeans(fpAnova,
#                      pairwise ~ condition|foreperiod,
#                      adjust = 'none')


fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')

# Compare variances
leveneTest(meanRT ~ condition * foreperiod,
           data = summaryData2)

#=========================== 2.2. FP x Accuracy by condition ===============================
acc_by_condition <- ggplot(data = summaryDataAll,
                           aes(x = foreperiod,
                               y = meanAcc,
                               color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 0.8, width = 0.1) +
  labs(x = "Foreperiod",
       y = "Mean accuracy",
       color = "condition",
       title = "Proportion of errors") +
  scale_color_manual(values = c("orange", "blue")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)))
ggsave("./Analysis/Plots/Acc_by_condition.png",
       acc_by_condition,
       width = 7.7,
       height = 5.8)

AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAll,
                   within = c("foreperiod", "condition"))

#================================= 3. Sequential effects ==============================
# Effects of previous orientation
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=prevOri)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = prevOri)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange'))

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=seqOri)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "Previous Orientation") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange'))


seqOriAnova <- aov_ez(id = "ID",
                   dv = "meanRT",
                   data = summaryData2,
                   within = c("foreperiod", "condition", "seqOri"))

#============================= 2.4. Learning effects ==============================
firstBlockData <- data2 %>%
  filter(block == '0', trial_bl %in% 1:80)

ggplot(data = firstBlockData,
       aes(x = trial_bl,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_wrap(~ foreperiod, nrow = 2, ncol = 1)


ggplot(data = firstBlockData,
       aes(x = action_trigger.rt,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean_cl_boot", geom = "errorbar", width = 0.2)

# Regression by block
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)

ggplot(data = data2,
       aes(x = trial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_wrap(~ foreperiod, nrow = 2, ncol = 1)

#==========================================================================================#
#====================================== 3. Mixed models ====================================
#==========================================================================================#

#=========================== 3.1. Foreperiod, condition and sequential effects =============================

# ============ 3.1.1. n-1 sequential effect =============
# 3.1.1.1 Fullest model
fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                 (1+numForeperiod*condition*numOneBackFP|ID),
               data = data,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'S',
               REML=TRUE)

fplmm1v2 <- lmer(RT ~ numForeperiod*condition*numOneBackFP + 
                   (1+numForeperiod*condition*numOneBackFP||ID),
                 data=data,
                 contrasts = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 REML=TRUE)

fplmm1v3 <- lmer(RT ~ numForeperiod*condition*numOneBackFP + 
                   (1+numForeperiod+condition+numOneBackFP||ID),
                 data=data,
                 contrasts = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 REML=TRUE)


# ============ 3.1.2. difference between the durations of FPn and FPn-1 as regressor ===========
# 3.1.1.1 Fullest model
fpDifflmm1 <- mixed(formula = RT ~ foreperiod * condition * oneBackFPDiff + 
                  (1+foreperiod*condition*oneBackFPDiff|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

fpDifflmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = data,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE)

#======================= 3.2. Foreperiod, condition and n-1 trial type =====================
# 3.2.1.1. Fullest model
oneBacktrialTypelmm1 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType + 
                                (1+foreperiod*condition*oneBacktrialType|ID),
                              data = data,
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'S',
                              REML=TRUE)

oneBacktrialTypelmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType + 
                                (1+foreperiod*condition*oneBacktrialType||ID),
                              data = data,
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'S',
                              REML=TRUE)

# 3.2.2
# 3.2.2.1. Remove correlations of mixed part
oneBacktrialTypeSeqlmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType * oneBackFP + 
                                   (1+foreperiod*condition*oneBacktrialType*oneBackFP||ID),
                                 data = data,
                                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                 progress = TRUE,
                                 expand_re = TRUE,
                                 method =  'S',
                                 REML=TRUE)

# 3.2.2.3. Remove interactions of mixed part
oneBacktrialTypeSeqlmm3 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType * oneBackFP + 
                                   (1+foreperiod+condition+oneBacktrialType+oneBackFP||ID),
                                 data = data,
                                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                 progress = TRUE,
                                 expand_re = TRUE,
                                 method =  'S',
                                 REML=TRUE)
#====================== 3.4. Foreperiod, counterbalance and block ========================
# 3.4.1. Fullest model
blockFplmm1 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                  (1+foreperiod*counterbalance*block|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

# 3.4.2. Remove correlations of mixed part
blockFplmm2 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                       (1+foreperiod*counterbalance*block||ID),
                     data = data,
                     control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                     progress = TRUE,
                     expand_re = TRUE,
                     method =  'S',
                     REML=TRUE)

# 3.4.3. Remove interactions of mixed part
blockFplmm3 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                       (1+foreperiod+counterbalance+block||ID),
                     data = data,
                     control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                     progress = TRUE,
                     expand_re = TRUE,
                     method =  'S',
                     REML=TRUE)
summary(blockFplmm3)$varcor
summary(blockFplmm3)
anova(blockFplmm3)


# 3.4.4. Plot
ggplot(data = data,
       aes(x = block,
           y = RT,
           color = counterbalance)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group = counterbalance)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("blue","orange")) +
  facet_wrap(~foreperiod)

# 3.5. Learning 





# Orientation
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color = orientation)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = orientation)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("deeppink3","chartreuse3"))

# Previous orientation
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color = prevOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = prevOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))

# Sequential effect by previous orientation


# Only action data
ggplot(data = filter(summaryData, condition == "action"),
       aes(x = foreperiod,
           y = meanRT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line" ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))

# Accuracy
ggplot(data = dataAll,
       aes(x = foreperiod,
           y = Acc)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group=1)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~ID)

# Fit by participant
ggplot(data = data,
       aes(x = foreperiod,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  facet_wrap(~ID) +
  scale_color_manual

# Fit by participant
ggplot(data = data,
       aes(x = foreperiod,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group=1)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  facet_wrap(~ID) +
  scale_color_manual(values = c("orange","blue"))
  #ylim(0.25,0.55)

# Boxplots
boxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanRT,
                       color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

boxplots


logBoxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanLogRT,
                       color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

logBoxplots


# Histograms
histograms <- ggplot(data=summaryData,
                     aes(x=meanRT,
                         color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)+
  theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

histograms <- ggplot(data=data,
                     aes(x=RT))+
  geom_histogram()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~ID)


ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanLogRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("orange", "blue"))


#================================== 4. Bayesian analysis ============================
bAnova <- anovaBF(meanRT ~ condition * foreperiod,
                  data = summaryData2)

bAnova/max(bAnova)


#====================== 5. Estimate n of trials to detect effect ======================
