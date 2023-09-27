

#================================================================================#
# Changes
# Manually set contrasts for anova
#================================================================================#

# Load necessary packages

# Read and process data
library(tidyverse)
library(broom)
library(magrittr)
library(data.table)

# Plotting
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)

# Linear models
library(car)
library(codingMatrices)
library(modelr)

# Mixed modeling
library(afex)
library(emmeans)
library(lme4)
library(performance)

# Bayesian analysis
library(BayesFactor)
library(bayestestR)

# Load data
source('./Analysis/Prepare_data_con.R')

# Create columns for FP length order and laterality
summaryData2 <- summaryData2 %>%
  arrange(ID, block) %>%
  group_by(ID) %>%
  mutate(fpOrder = ifelse(foreperiod[1] == 1000, "short-long", "long-short")) %>%
  ungroup()
  # ungroup() %>%
  # mutate(laterality = rep(c("right", "left", "right", "right", "right", "right", NaN, "right", "right", "left", "right",
  #                       "right", "right", "right", "right", "right", "right", "right", "right", "right",
  #                       NaN, "right", "right", "right", "right", "left", "right", NaN, "right"), each = 4))

# Check for influence of laterality
# summaryData2 <- summaryData2 %>%
#   filter(laterality != "left")

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
ggplot(data=summaryData2,
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

# By FP order
ggplot(data = summaryData2,
       aes(x = block,
           y = meanRT,
           color = fpOrder)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = fpOrder)) +
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('blue','orange'))

# Distribution of data
dataHists <- ggplot(data=summaryData2,
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
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8))) +
  labs(x = "External fixation duration",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/extfixduration.png",
       width = 13.4,
       height = 10)

# Check for influence of latency of action key press on RT
ggplot(data=filter(data,condition=='action',action_trigger.rt < 5.0),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8))) +
  labs(x = "Action trigger delay",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/actiontrigpress.png",
       width = 13.4,
       height = 10)


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

jpeg("./Analysis/Plots/BF_seq.jpeg", width = 1200, height = 1000)
plot(srange, fp_bfs,
     xlab = "n of participants",
     ylab = "Bayes Factor",
     main = "Bayes Factors for FP effect compared between conditions",
     cex.lab = 2,
     cex.axis = 2,
     cex.main = 2.5)
lines(srange, fp_bfs)
dev.off()

# % of extreme RTs removed
(ntrials_before_extrem - ntrials_after_extrem)/ntrials_before_extrem * 100

# % of RTs removed due to trimming by participant
n_notrim <- data %>%
  group_by(ID) %>%
  summarise(n_trials_no_trim = n())

n_trim <- data2 %>%
  group_by(ID) %>%
  summarise(n_trials_trim = n())

n_trials_comp <- left_join(n_notrim, n_trim) %>%
  mutate(trimmed = n_trials_no_trim - n_trials_trim,
         percent_trimmed = (trimmed/n_trials_no_trim) * 100) %>%
  summarise(meanT = mean(percent_trimmed), sdT = sd(percent_trimmed))

#==========================================================================================#
#==================================== 1. Descriptives ======================================
#==========================================================================================#

#======================================== Plots ============================================
# main effect of condition
ggplot(data = summaryData2,
       aes(x = condition,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = 'mean', geom = 'point') +
  stat_summary(fun = 'mean', geom = 'line' ) +
  stat_summary(fun.data = 'mean_cl_boot', width = 0.2, geom = 'errorbar') +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))

# RT by FP and condition
RT_by_condition_part <- ggplot(data = summaryData2,
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
                RT_by_condition_part)

RT_by_condition <- ggplot(data = summaryData2 %>%
                            group_by(ID, foreperiod, condition) %>%
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", linewidth = 1.2, width = 0.1, geom = "errorbar") +
  labs(title = "RT",
       x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt")) +
  scale_color_manual(values = c("orange","blue"),
                     label = c("External", "Action"))
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.pdf",
                RT_by_condition,
                width = 7.7,
                height = 5.8)

acc_by_condition <- ggplot(data = summaryDataAcc %>%
                             group_by(ID, foreperiod, condition) %>%
                             summarise(meanAcc = mean(meanAcc)),
                           aes(x = foreperiod,
                               y = meanAcc,
                               color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", linewidth = 1.2, width = 0.1, geom = "errorbar") +
  labs(x = "FP (s)",
       y = "Mean prop. correct",
       color = "Condition",
       title = "Accuracy") +
  scale_color_manual(values = c("orange", "blue"),
                     label = c("External", "Action")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"))
ggsave("./Analysis/Plots/Acc_by_condition.pdf",
       acc_by_condition,
       width = 7.7,
       height = 5.8)

error_by_condition <- ggplot(data = summaryDataAcc %>%
                             group_by(ID, foreperiod, condition) %>%
                             summarise(errorRate = mean(errorRate)),
                           aes(x = foreperiod,
                               y = errorRate,
                               color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", linewidth = 1.2, width = 0.1, geom = "errorbar") +
  labs(x = "FP (s)",
       y = "Mean Error Rate",
       color = "Condition",
       title = "Error Rate") +
  scale_color_manual(values = c("orange", "blue"),
                     label = c("External", "Action")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt"))
ggsave("./Analysis/Plots/error_by_condition.pdf",
       error_by_condition,
       width = 7.7,
       height = 5.8)

# RT and accuracy in single panel
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(acc_by_condition + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")


# Visualize
grid.arrange(RT_by_condition + theme(legend.position = "none",
                                     axis.text = element_text(size = rel(1.2)),
                                     axis.title = element_text(size = rel(1.4)),
                                     plot.title = element_text(size = rel(1.5))),
             acc_by_condition + theme(legend.position = "none",
                                      axis.text = element_text(size = rel(1.2)),
                                      axis.title = element_text(size = rel(1.4)),
                                      plot.title = element_text(size = rel(1.5))),
             cond_legend,
             nrow = 1,
             widths = c(4/9, 4/9, 1/9))

# Save plots
comb_plots <- arrangeGrob(RT_by_condition + theme(legend.position = "none",
                                                  axis.text = element_text(size = rel(1.2)),
                                                  axis.title = element_text(size = rel(1.4)),
                                                  plot.title = element_text(size = rel(1.5))),
                          acc_by_condition + theme(legend.position = "none",
                                                   axis.text = element_text(size = rel(1.2)),
                                                   axis.title = element_text(size = rel(1.4)),
                                                   plot.title = element_text(size = rel(1.5))),
                          cond_legend,
                          nrow = 1,
                          widths = c(4/9, 4/9, 1/9))

ggsave("./Analysis/Plots/comb_plots.pdf",
       comb_plots,
       width = 20,
       height = 11.11,
       unit = "cm")

# RT and error rate in single panel
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(error_by_condition + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")

xaxis_title <- text_grob(error_by_condition$labels$x,
                         just = "top",
                         size = (error_by_condition + 
                                   theme(axis.title = element_text(size = rel(1.4))))$theme$axis.title$size * 11) # 11 is the base size in theme_grey


xaxis_title_margin <- unit(2, "line")


# Visualize
grid.arrange(arrangeGrob(RT_by_condition + theme(legend.position = "none",
                                     axis.text = element_text(size = rel(1.2)),
                                     axis.title = element_text(size = rel(1.4)),
                                     plot.title = element_text(size = rel(1.5)),
                                     axis.title.x = element_blank()),
             error_by_condition + theme(legend.position = "none",
                                      axis.text = element_text(size = rel(1.2)),
                                      axis.title = element_text(size = rel(1.4)),
                                      plot.title = element_text(size = rel(1.5)),
                                      axis.title.x = element_blank()),
             cond_legend,
             nrow = 1,
             widths = c(4/9, 4/9, 1/9)),
             xaxis_title,
             heights = unit.c(unit(1, "null"),
                              grobHeight(xaxis_title) + xaxis_title_margin),
             nrow = 2)

# Save plots
rt_error_plots <- arrangeGrob(arrangeGrob(RT_by_condition + theme(legend.position = "none",
                                                                  axis.text = element_text(size = rel(1.2)),
                                                                  axis.title = element_text(size = rel(1.4)),
                                                                  plot.title = element_text(size = rel(1.5)),
                                                                  axis.title.x = element_blank()),
                                          error_by_condition + theme(legend.position = "none",
                                                                     axis.text = element_text(size = rel(1.2)),
                                                                     axis.title = element_text(size = rel(1.4)),
                                                                     plot.title = element_text(size = rel(1.5)),
                                                                     axis.title.x = element_blank()),
                                          cond_legend,
                                          nrow = 1,
                                          widths = c(4/9, 4/9, 1/9)),
                              xaxis_title,
                              heights = unit.c(unit(1, "null"),
                                               grobHeight(xaxis_title) + xaxis_title_margin),
                              nrow = 2)

ggsave("./Analysis/Plots/rt_error_plots.pdf",
       rt_error_plots,
       width = 20,
       height = 11.11,
       unit = "cm")

#==================================== 1.2. Tables ================================
# By foreperiod and condition
meandata <- summaryData2 %>%
  group_by(condition, foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandata  

# By foreperiod
meandataFP <- summaryData2 %>%
  group_by(foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandataFP  


# By condition
#==========================================================================================#
#==================================== 2. Basic models ======================================
#==========================================================================================#

# Set constrasts for variables used in ANOVAs
contrasts(summaryData$foreperiod) <- c(-1/2, 1/2)
contrasts(summaryData$condition) <- c(-1/2, 1/2)

contrasts(summaryData2$foreperiod) <- c(-1/2, 1/2)
contrasts(summaryData2$condition) <- c(-1/2, 1/2)

#==================== 2.1. FP x RT by condition ======================
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


# Bootstrap
library(marginaleffects)


fpComp <- avg_comparisons(fpAnova$lm, variables = c("condition", "foreperiod"))

fpComp <- avg_comparisons(fpAnova, by = "foreperiod", variables = "condition")

fpComp |> inferences(method = "boot")


# Try transformations - invRT
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

# Try transformations - logRT
logfpAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"),
                     anova_table = list(es = "pes"))

### Check assumptions

# Sphericity
check_sphericity(logfpAnova)

# Normality of residuals
is_norm <- check_normality(logfpAnova)

plot(is_norm)

plot(is_norm, type = 'qq')

plot(is_norm, type = 'qq', detrend = TRUE)

# The log-transform does not solve the problem either



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

# Anova
varAnova <- aov_ez(data = summaryData2,
                   id = "ID",
                   dv = "varRT",
                   within = c("foreperiod", "condition"),
                   anova_table = list(es = "pes"))

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = sqrt(varRT),
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")

#========================== 2.2. Orientation sequential effects ============================
# Orientation
ggplot(data = summaryData2,
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

# Orientation and condition
ggplot(data = summaryData2,
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
        axis.title = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5))) +
  labs(color = "Orientation",
       x = "Foreperiod",
       y = "mean RT") +
  facet_wrap(~ condition) +
  scale_color_manual(values = c("deeppink3","chartreuse3"))
ggsave("./Analysis/Plots/RT_orientation_condition.png",
       width = 13.4,
       height = 10)

# Previous orientation
ggplot(data = summaryData2,
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
  scale_color_manual(values = c("lightgoldenrod1","indianred2")) +
  facet_wrap(~ condition)

# Repetition/alternation
ggplot(data = summaryData2 %>%
         filter(!is.na(seqOri)),
       aes(x = foreperiod,
           y = meanRT,
           color = seqOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2")) +
  facet_wrap(~ condition)


# Repetition/alternation
ggplot(data = summaryData2 %>%
         filter(!is.na(seqOri)),
       aes(x = foreperiod,
           y = meanRT,
           color = seqOri)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = seqOri)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.5))) +
  labs(color = "Previous orientation",
       x = "Foreperiod",
       y = "mean RT") +
  facet_wrap(~ condition) +
  scale_color_manual(values = c("lightgoldenrod1","indianred2"))

ggsave("./Analysis/Plots/RT_seqalt_condition.png",
       width = 13.4,
       height = 10)

#=========================== 2.3. FP x Accuracy by condition ===============================


AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAcc,
                   within = c("foreperiod", "condition"),
                   anova_table = list(es = "pes"))


check_sphericity(AccAnova)

is_norm_acc <- check_normality(AccAnova)

# Visualize
groupedAcc <- summaryDataAcc %>%
  group_by(ID, foreperiod, condition) %>%
  summarise(meanAcc = mean(meanAcc))


hist(groupedAcc$meanAcc)

plot(is_norm_acc)

plot(is_norm_acc, type = "qq")

plot(is_norm_acc, type = "qq", detrend = TRUE)

friedman.test(y = groupedAcc$meanAcc, groups = c(groupedAcc$condition, groupedAcc$foreperiod), blocks = groupedAcc$ID)


library(permuco)

aovperm(meanAcc ~ foreperiod * condition + Error(ID/(foreperiod * condition)), data = summaryDataAcc)

#============================= 2.4. Learning effects ==============================
firstBlockData <- data2 %>%
  filter(block == '0', trial_bl %in% 1:80)


ggplot(data = firstBlockData,
       aes(x = trial_bl,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c('orange', 'blue'))# +
  #facet_wrap(~ foreperiod, nrow = 2, ncol = 1)


ggplot(data = firstBlockData,
       aes(x = action_trigger.rt,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean_cl_boot", geom = "errorbar", width = 0.2)

# Binned trials
data2 <- data2 %>%
  group_by(ID) %>%
  mutate(binTrial = ntile(trial, n = 16)) %>%
  ungroup()


ggplot(data = data2,
       aes(x = binTrial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  scale_color_manual(values = c('orange', 'blue'))

firstBlockData <- firstBlockData %>%
  group_by(ID) %>%
  mutate(binTrial = ntile(trial, n = 8)) %>%
  ungroup()

ggplot(data = firstBlockData,
       aes(x = binTrial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  scale_color_manual(values = c('orange', 'blue'))


# Regression by block
blocklm <- lm(meanRT ~ foreperiod * condition * block,
              data = summaryData2)

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



#================================== 3. Bayesian analysis ============================
bAnova <- anovaBF(meanRT ~ condition * foreperiod,
                  data = summaryData2)

bAnova/max(bAnova)

