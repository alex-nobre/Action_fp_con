

# Load packages

# Data processing and plotting
library(magrittr)
library(tidyverse)
library(lattice)
library(gridExtra)
library(data.table)

# Simple models
library(car)
library(janitor)

# Mixed-effects modeling
library(afex)
library(emmeans)
library(lme4)
library(MuMIn)
library(buildmer)
library(broom.mixed)

# Bayesian models
library(brms)
library(bayestestR)
library(BayesFactor)

# Assess models and results
library(effects)
library(ggeffects)
library(performance)
library(knitr)
library(kableExtra)
library(sjPlot)
library(prediction)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 


#================================== 0. Read data ================================
# Create dataset
source('./Analysis/Prepare_data_con.R')

# Set contrasts
contrasts(data$foreperiod) <- c(-1/2, 1/2)
contrasts(data$logFP) <- c(-1/2, 1/2)
contrasts(data$condition) <- c(-1/2, 1/2)
contrasts(data$prevOri) <- c(-1/2, 1/2)
contrasts(data$seqOri) <- c(-1/2, 1/2)

contrasts(data2$foreperiod) <- c(-1/2, 1/2)
contrasts(data2$logFP) <- c(-1/2, 1/2)
contrasts(data2$condition) <- c(-1/2, 1/2)
contrasts(data2$prevOri) <- c(-1/2, 1/2)
contrasts(data2$seqOri) <- c(-1/2, 1/2)

#==========================================================================================#
#================================= 1. Explore individual data ==============================
#==========================================================================================#

# Functions
hist_resid <- function(M,ptitle='Residuals') {
  d <- data.frame(resid=residuals(M)) 
  d  %>% ggplot(aes(x=resid)) + 
    geom_histogram(aes(y=after_stat(density)), bins=75, color='black', fill='grey') + 
    geom_density(color='darkred') + 
    ggtitle(ptitle) -> pl
  return(pl)
}

fitstats = function(M,mname='M') {
  QQ<-qqnorm(residuals(M), plot.it=FALSE)
  R2qq <- cor(QQ$x,QQ$y)^2
  dfqq = data.frame(stat='R2qq', V1=R2qq)
  r2tab <- r.squaredGLMM(M)  %>% 
    t  %>% as.data.frame  %>% rownames_to_column(var='stat')  %>% 
    rbind(.,dfqq)
  r2tab$stat = c("$R^2_m$","$R^2_c$",'$R^2_{qq}$' )
  colnames(r2tab) <- c('stat',mname)
  return(r2tab)
}

# Plot RT by FP by participant and model using complete pooling 
FPfitAll=lm(meanRT ~ foreperiod,
            data=summaryData)

fit.params=tidy(FPfitAll)

summary(FPfitAll)


ggplot(data=summaryData,
       aes(x=foreperiod,
           y=meanRT)) +
  stat_summary(fun="mean", geom="point", size=1.5)+
  geom_abline(intercept=fit.params$estimate[1],
              slope=fit.params$estimate[2],
              color="blue")+
  facet_wrap(~ ID, ncol=6)


# Plot RT by FP by participant and model using individual data (no pooling)
dataGroupedByRT <- summaryData %>% 
  group_by(ID,foreperiod) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

data.no_pooling <- dataGroupedByRT %>%
  select(-foreperiod) %>%
  group_by(ID) %>%
  nest(data = c(numForeperiod, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numForeperiod, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, term, estimate) %>%
  complete(ID, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(ID=id,
         numForeperiod=num_foreperiod)


ggplot(data = dataGroupedByRT,
       aes(x = numForeperiod, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numForeperiod),
              color = "blue") +
  geom_point() +
  facet_wrap(~ID, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))

fp_no_pooling <- data.no_pooling$numForeperiod

# Compare results to see how much they differ
data_grouped_by_fp <- data %>%
  group_by(ID, foreperiod) %>%
  summarise(meanRT=mean(RT)) %>%
  ungroup()

fit_partial_pooling <- lmer(formula = RT ~ foreperiod + 
                              (1 + foreperiod|ID),
                            data = data)

data_partial_pooling <- fit_partial_pooling %>%
  augment() %>%
  select(ID, foreperiod, RT, .fitted) %>%
  rename(fitted=.fitted)


#==========================================================================================#
#====================================== 2. Prepare model ====================================
#==========================================================================================#

#=========================== 2.1. Choose dependent variable =============================
# To choose the data transformation that leads to the optimal random effects structure, we fit models including only 
# random intercepts and compare R2 and residuals

# Fit models with RT and inverse RT without trimming


fplmm <- mixed(formula = RT ~ foreperiod * condition + 
                  (1|ID),
                data = data,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

# Now we run the same model with inverse RT and logRT as outcomes
invfplmm <- mixed(formula = invRT ~ foreperiod * condition+ 
                     (1|ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

logfplmm <- mixed(formula = logRT ~ foreperiod * condition + 
                     (1|ID),
                   data = data,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")


# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm), data$logRT)^2)

# invRt explains the most variance


# Check normality of residuals
par(mfrow=c(1,3))
qqnorm(resid(fplmm),
       main="Normal q-qplot fplmm")
qqnorm(resid(invfplmm),
       main="Normal q-qplot invfplmm")
qqnorm(resid(logfplmm),
       main="Normal q-qplot logfplmm")

par(graphical_defaults)

# All models show moderate departures from normality, with logRT outperforming the others

# Plot residuals
par(mfrow=c(3,1))
plot(resid(fplmm), fitted(fplmm),
     main="Residuals fplmm")
plot(resid(invfplmm), fitted(invfplmm),
     main="Residuals invfplmm")
plot(resid(logfplmm), fitted(logfplmm),
     main="Residuals logfplmm")

par(graphical_defaults)
# There are no clear correlations

# Residual histograms
grid.arrange(hist_resid(fplmm, 'RT'),
             hist_resid(invfplmm, '1/RT'),
             hist_resid(logfplmm, 'logRT'),
             ncol=1)

# All appear to be relatively normally distributed, with the smallest skew for logRT

# Fit models with RT, inverse RT and log RT with trimming
trimfplmm <- mixed(formula = RT ~ foreperiod * condition + 
                      (1|ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")


# Now we run the same model with inverse RT and logRT as outcomes
triminvfplmm <- mixed(formula = invRT ~ foreperiod * condition + 
                         (1|ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")


trimlogfplmm <- mixed(formula = logRT ~ foreperiod * condition + 
                         (1|ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")


# Amount of variance accounted for by the model
var <- rbind(var,
             data.frame(dataset = 'trim',
                        'RT' = cor(fitted(trimfplmm), data2$RT)^2,
                        'invRT' = cor(fitted(triminvfplmm), data2$invRT)^2,
                        'logRT'= cor(fitted(trimlogfplmm), data2$logRT)^2))


var
# The invRT model accounts for the most variance. Trimming appears to make a small
# difference

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm),
       main="Normal q-qplot fplmm")
qqnorm(resid(invfplmm),
       main="Normal q-qplot invfplmm")
qqnorm(resid(logfplmm),
       main="Normal q-qplot logfplmm")
qqnorm(resid(trimfplmm),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality, and trimming does seem to help a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm), fitted(fplmm),
     main="Residuals RT")
plot(resid(trimfplmm), fitted(trimfplmm),
     main="Residuals trimmed RT")
plot(resid(invfplmm), fitted(invfplmm),
     main="Residuals 1/RT")
plot(resid(triminvfplmm), fitted(triminvfplmm),
     main="Residuals trimmed 1/RT")
plot(resid(logfplmm), fitted(logfplmm),
     main="Residuals logRT")
plot(resid(trimlogfplmm), fitted(trimlogfplmm),
     main="Residuals trimmed logRT")
par(graphical_defaults)

# LIttle difference after trimming, although outliers are rarer; logRT still appears to perform better

grid.arrange(hist_resid(fplmm, 'RT'),
             hist_resid(invfplmm, '1/RT'),
             hist_resid(logfplmm, 'logRT'),
             hist_resid(trimfplmm, 'trimmed RT'),
             hist_resid(triminvfplmm, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)

# logRT outperforms the others, and trimming seems to help
R2table <- fitstats(fplmm, 'RT') %>%
  plyr::join(., fitstats(invfplmm, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm, 'trim log(RT)'), by='stat') %>%
  knitr::kable(digits=4, format = "pipe")

# Indeed, the model using trimmed invRT has the highest R2 values, although it is a small difference.
# All R2 values are very low 

# Variance explained is higher with trimming, although this is probably not significant
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT explains a larger part of the variance, but does less well on residuals

# Since the differences in R2 are small, but the differences in residuals appear considerable, and due to 
# ease of interpretation, we use logRT

#================================= 2.2. Find random effects structure ==========================

trimlog <- buildmer(logRT ~ foreperiod * condition + 
                            (1 + foreperiod * condition |ID), 
                          data=data2,
                          buildmerControl = list(family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))

isSingular(trimlog)
formula(trimlog)

# Systematic comparisons betweeen lmm's via BIC

# Model obtained with buildmer
trimlog <- mixed(logRT ~ 1 + foreperiod + condition + foreperiod:condition + 
                    (1 + condition + foreperiod + foreperiod:condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

anova(trimlog)

# Random-intercept only model
trimlog1v2 <- mixed(logRT ~ 1 + foreperiod * condition + 
                         (1 | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)


# Compare BICs and AICs
BIC(trimlog1, trimlog1v2)

AIC(trimlog1, trimlog1v2)

cor(fitted(trimlog1), data2$logRT)^2
cor(fitted(trimlog1v2), data2$logRT)^2

# Full model explains a larger part of the variance

#==============================================================================================#
#==================================== 3. Model assessment ======================================
#==============================================================================================#

#==================================== 3.1. Check effects =======================================
trimlog <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + 
                   (1 + condition + foreperiod + foreperiod:condition | ID),
                 data = data2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod")


isSingular(trimlog)

anova(trimlog)

#================== 3.3. Compare dependent variables using random-effects structure ==================
trimfplmm1 <- mixed(formula = RT ~ condition * foreperiod
                      (1 + condition | ID),
                    data = data2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")


triminvfplmm1 <- mixed(formula = invRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

trimlogfplmm1 <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data = data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# Compare model R2 and residuals
# Amount of variance accounted for by the model
var <- data.frame(dataset = 'no trim',
                  'RT' = cor(fitted(fplmm1), data$RT)^2,
                  'invRT' = cor(fitted(invfplmm1), data$invRT)^2,
                  'logRT' = cor(fitted(logfplmm1), data$logRT)^2)

var <- rbind(var,
             data.frame(dataset = 'trim',
                        'RT' = cor(fitted(trimfplmm1), data2$RT)^2,
                        'invRT' = cor(fitted(triminvfplmm1), data2$invRT)^2,
                        'logRT'= cor(fitted(trimlogfplmm1), data2$logRT)^2))

var
# Again, the log model accounts (barely) for a larger amount of the variance

# check normality
par(mfrow=c(2,3))
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm")
par(graphical_defaults)

# The log plots show the least deviations from normality, and trimming does seem to help a bit

# Plot residuals
par(mfrow=c(3,2))
plot(resid(fplmm1), fitted(fplmm1),
     main="Residuals RT")
plot(resid(invfplmm1), fitted(invfplmm1),
     main="Residuals 1/RT")
plot(resid(logfplmm1), fitted(logfplmm1),
     main="Residuals logRT")
plot(resid(trimfplmm1), fitted(trimfplmm1),
     main="Residuals trimmed RT")
plot(resid(triminvfplmm1), fitted(triminvfplmm1),
     main="Residuals trimmed 1/RT")
plot(resid(trimlogfplmm1), fitted(trimlogfplmm1),
     main="Residuals trimmed logRT")
par(graphical_defaults)

grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)


R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)


# By most measures, logRT with trimming still yields the best fit

# Compare lm and lmm models
trimlogfplm <- lm(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                    numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP,
                  data = data2)

trimlogfplmmML <- mixed(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data=data2,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        REML=FALSE,
                        return = "merMod",
                        check_contrasts = FALSE)

trimlogfplmmML <- lmer(logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                         numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                         (1 + condition | ID),
                       data=data2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       method =  'KR',
                       REML=FALSE)

anova(trimlogfplmmML, trimlogfplm)

#============ 3.5. Hierarchical entry ===============
h_trimlogfplmm1 <- mixed(logRT ~ 1 + numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method = 'KR',
                         REML=TRUE,
                         return = "merMod")

h_trimlogfplmm2 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm3 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm4 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm5 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm6 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

h_trimlogfplmm7 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)


h_trimlogfplmm8 <- mixed(logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + condition + 
                           condition:numForeperiod + condition:numOneBackFP + condition:numOneBackFP:numForeperiod + (1 + condition | ID),
                         data=data2,
                         control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                         progress = TRUE,
                         expand_re = TRUE,
                         method =  'KR',
                         REML=TRUE,
                         return = "merMod",
                         check_contrasts = FALSE)

anova(h_trimlogfplmm1, h_trimlogfplmm2, h_trimlogfplmm3, h_trimlogfplmm4, 
      h_trimlogfplmm5, h_trimlogfplmm6, h_trimlogfplmm7, h_trimlogfplmm8)


#==================== 3.3. Run model separately for action and external conditions ===================

#======== 3.3.1. External ============#

# FP as categorical
trimlogfplmmext2 <- mixed(logRT ~ 1 + foreperiod + 
                            numOneBackFP + foreperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='external',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmmext2)


#======== 3.3.2. Action ============#
trimlogfplmm1act <- mixed(logRT ~ 1 + numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP +
                            (1 | ID),
                          data=data2[data2$condition=='action',],
                          control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                          progress = TRUE,
                          expand_re = TRUE,
                          method =  'KR',
                          REML=TRUE,
                          return = "merMod",
                          check_contrasts = FALSE)

summary(trimlogfplmm1act)


#==============================================================================================#
#================================== 4. Choose distribution ====================================
#==============================================================================================#

#======================== 4.1. Visualize distributions for each variable =======================
# Pooled
RTHistograms <- ggplot(data=data2,
                       aes(x=RT))+
  geom_histogram()
RTHistograms

invRTHistograms <- ggplot(data=data2,
                          aes(x=invRT)) +
  geom_histogram()
invRTHistograms

logRTHistograms <- ggplot(data=data2,
                          aes(x=logRT)) +
  geom_histogram()
logRTHistograms


# By participant
indRTHistograms <- ggplot(data=data2,
                          aes(x=RT))+
  geom_histogram()+
  facet_wrap(~ID)
indRTHistograms

indinvRTHistograms <- ggplot(data=data2,
                             aes(x=invRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indinvRTHistograms

indlogRTHistograms <- ggplot(data=data2,
                             aes(x=logRT)) +
  geom_histogram() +
  facet_wrap(~ID)
indlogRTHistograms

# Try models
trimlogfpgauss <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data = data2,
                        #family=gaussian(link = "identity"),
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimlogfpgauss)

trimlogfpinvgauss <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                             numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                             (1 + condition | ID),
                           data = data2,
                           family=inverse.gaussian(link = "identity"),
                           control = lmerControl(optimizer = c("nloptwrap"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                           progress = TRUE,
                           expand_re = TRUE,
                           method =  'KR',
                           return = "merMod")

summary(trimlogfpinvgauss)

trimlogfpgamma <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          numOneBackFP + numForeperiod:numOneBackFP + condition:numOneBackFP + 
                          (1 + condition | ID),
                        data = data2,
                        family=Gamma(link = "identity"),
                        control = lmerControl(optimizer = c("nloptwrap"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        return = "merMod")

summary(trimlogfpgamma)

# Compare visualizations
ggpredict(model = trimlogfpgauss,
          terms = "numForeperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpinvgauss,
          terms = "numForeperiod",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgamma,
          terms = 'numForeperiod',
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpinvgauss,
          terms = "condition",
          type = 'fe') %>%
  plot()

ggpredict(model = trimlogfpgamma,
          terms = 'condition',
          type = 'fe') %>%
  plot()


# Compare performance across models
trimlogfpgauss %>%
  check_model()

trimlogfpinvgauss %>%
  check_model()

trimlogfpgamma %>%
  check_model()

# The gaussian model seems to provide a better fit than the alternatives

#==============================================================================================#
#============================== 5. Bayesian mixed models using brms ============================
#==============================================================================================#

#========================== 5.1. For Log RT ==============================
#======= 5.1.1 Models using effects from previous experiment as priors

# Baseado no exp 0 v1
# prior1 <- c(set_prior('normal(2.55,0.009)', class = 'Intercept'),
#             set_prior("normal(-0.03, 0.002)", class = "b", coef = "foreperiod2"),
#             set_prior("normal(-0.02, 0.002)", class = "b", coef = "foreperiod3"),
#             set_prior("normal(0.01, 0.002)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(0.02, 0.002)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(0.01, 0.007)", class = "b", coef = "condition1"),
#             set_prior("normal(-0.02, 0.005)", class = "b", coef = "foreperiod2:oneBackFP2"),
#             set_prior("normal(-0.03, 0.005)", class = "b", coef = "foreperiod3:oneBackFP2"),
#             set_prior("normal(-0.03, 0.005)", class = "b", coef = "foreperiod2:oneBackFP3"),
#             set_prior("normal(-0.04, 0.005)", class = "b", coef = "foreperiod3:oneBackFP3"),
#             set_prior("normal(-0.01, 0.004)", class = "b", coef = "foreperiod2:condition1"),
#             set_prior("normal(-0.05, 0.004)", class = "b", coef = "foreperiod3:condition1"),
#             set_prior("normal(0.003, 0.004)", class = "b", coef = "oneBackFP2:condition1"),
#             set_prior("normal(0.002, 0.004)", class = "b", coef = "oneBackFP3:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod2:oneBackFP2:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod3:oneBackFP2:condition1"),
#             set_prior("normal(-0.02, 0.009)", class = "b", coef = "foreperiod2:oneBackFP3:condition1"),
#             set_prior("normal(-0.03, 0.009)", class = "b", coef = "foreperiod3:oneBackFP3:condition1")
# )
prior1 <- c(set_prior('normal(2.55,0.009)', class = 'Intercept'),
            set_prior("normal(-0.03, 0.002)", class = "b", coef = "foreperiod2"),
            set_prior("normal(-0.02, 0.002)", class = "b", coef = "foreperiod3"),
            set_prior("normal(0.01, 0.002)", class = "b", coef = "oneBackFP2"),
            set_prior("normal(0.02, 0.002)", class = "b", coef = "oneBackFP3"),
            set_prior("normal(0.01, 0.007)", class = "b", coef = "condition1")
)

# Gaussian distr
b_one_back_fp <- brm(formula = logRT ~ foreperiod * condition * oneBackFP + 
                       (1+numForeperiod*condition*numOneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior1,
                     control = list(adapt_delta = 0.9),
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 12000,
                     cores = -1)

saveRDS(b_one_back_fp, "b_one_back_fp.rds")

#b_one_back_fp <- readRDS('./Analysis/b_one_back_fp.rds')


ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_point(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')



# Inverse gaussian distr
b_one_back_fp_invgaus <- brm(formula = RT ~ numForeperiod * condition * numOneBackFP + 
                               (1+numForeperiod*condition*numOneBackFP|ID),
                             data = data2,
                             family = gaussian(),
                             prior = prior1,
                             save_all_pars = TRUE,
                             warmup = 2000,
                             iter = 10000)


b_one_back_fp_null <- brm(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                            numOneBackFP + numForeperiod:numOneBackFP + 
                            (1 + condition + numForeperiod + numOneBackFP + 
                               condition:numForeperiod + numForeperiod:numOneBackFP| ID),
                          data = data2,
                          family = gaussian(),
                          prior = prior1,
                          save_all_pars = TRUE,
                          control = list(adapt_delta = 0.9),
                          warmup = 2000,
                          iter = 10000,
                          file="b_one_back_fp_null.rds")



bf_one_back_brm <- bayes_factor(b_one_back_fp, b_one_back_fp_null)


#========== 5.1.2 Using a vague prior
b_one_back_fp_vagueprior <- brm(formula = logRT ~ numForeperiod * condition * numOneBackFP + 
                                  (1+numForeperiod*condition*numOneBackFP|ID),
                                data = data2,
                                family = gaussian(),
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.9),
                                warmup = 2000,
                                iter = 12000,
                                cores = -1)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)

#===================================== 5.2. For linear RT ======================================
#======= 5.2.1 Models using effects from previous experiment as priors

# Baseado no exp 0 v1
# prior2 <- c(set_prior('normal(368.41,8.34)', class = 'Intercept'),
#             set_prior("normal(-27.82, 1.78)", class = "b", coef = "foreperiod2"),
#             set_prior("normal(-24.96, 1.78)", class = "b", coef = "foreperiod3"),
#             set_prior("normal(8.57, 1.77)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(17.22, 1.77)", class = "b", coef = "oneBackFP2"),
#             set_prior("normal(12.6, 6.27)", class = "b", coef = "condition1"),
#             set_prior("normal(-14.39, 4.38)", class = "b", coef = "foreperiod2:oneBackFP2"),
#             set_prior("normal(-24.31, 4.37)", class = "b", coef = "foreperiod3:oneBackFP2"),
#             set_prior("normal(-26.36, 4.34)", class = "b", coef = "foreperiod2:oneBackFP3"),
#             set_prior("normal(-37.74, 4.37)", class = "b", coef = "foreperiod3:oneBackFP3"),
#             set_prior("normal(-15.27, 3.55)", class = "b", coef = "foreperiod2:condition1"),
#             set_prior("normal(-22.61, 3.56)", class = "b", coef = "foreperiod3:condition1"),
#             set_prior("normal(2.92, 3.55)", class = "b", coef = "oneBackFP2:condition1"),
#             set_prior("normal(1.45, 3.54)", class = "b", coef = "oneBackFP3:condition1"),
#             set_prior("normal(-12.49, 8.75)", class = "b", coef = "foreperiod2:oneBackFP2:condition1"),
#             set_prior("normal(-21.64, 8.74)", class = "b", coef = "foreperiod3:oneBackFP2:condition1"),
#             set_prior("normal(-20.01, 8.70)", class = "b", coef = "foreperiod2:oneBackFP3:condition1"),
#             set_prior("normal(-23.69, 8.75)", class = "b", coef = "foreperiod3:oneBackFP3:condition1")
# )
prior2 <- c(set_prior('normal(368.41,8.34)', class = 'Intercept'),
            set_prior("normal(-27.82, 1.78)", class = "b", coef = "foreperiod2"),
            set_prior("normal(-24.96, 1.78)", class = "b", coef = "foreperiod3"),
            set_prior("normal(8.57, 1.77)", class = "b", coef = "oneBackFP2"),
            set_prior("normal(17.22, 1.77)", class = "b", coef = "oneBackFP3"),
            set_prior("normal(12.6, 6.27)", class = "b", coef = "condition1")
)


# Gaussian distr
b_one_back_fp <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                       (1+foreperiod*condition*oneBackFP|ID),
                     data = data2,
                     family = gaussian(),
                     prior = prior2,
                     control = list(adapt_delta = 0.9),
                     save_all_pars = TRUE,
                     warmup = 2000,
                     iter = 12000,
                     cores = -1)

saveRDS(b_one_back_fp, "b_one_back_fp.rds")

#b_one_back_fp <- readRDS('./Analysis/b_one_back_fp.rds')


ggplot() +
  stat_summary(fun='mean', geom='point', data=data2, aes(x=numForeperiod, y=logRT)) +
  geom_point(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_line(data=fp_effect_df, aes(x=numForeperiod, y=fit), color='red') +
  geom_ribbon(data=fp_effect_df, aes(x=numForeperiod, ymin=lower, ymax=upper), alpha=0.3) +
  labs(x='Foreperiod (continuous)', y='logRT')



# Inverse gaussian distr
b_one_back_fp_invgaus <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                               (1+foreperiod*condition*oneBackFP|ID),
                             data = data2,
                             family = gaussian(),
                             prior = prior1,
                             save_all_pars = TRUE,
                             warmup = 2000,
                             iter = 10000)


b_one_back_fp_null <- brm(formula = RT ~ 1 + condition + foreperiod + condition:foreperiod + 
                            oneBackFP + foreperiod:oneBackFP + 
                            (1 + condition + foreperiod + oneBackFP + 
                               condition:numForeperiod + numForeperiod:numOneBackFP| ID),
                          data = data2,
                          family = gaussian(),
                          prior = prior1,
                          save_all_pars = TRUE,
                          control = list(adapt_delta = 0.9),
                          warmup = 2000,
                          iter = 10000,
                          file="b_one_back_fp_null.rds")



bf_one_back_brm <- bayes_factor(b_one_back_fp, b_one_back_fp_null)


#========== 5.1.2 Using a vague prior
b_one_back_fp_vagueprior <- brm(formula = RT ~ foreperiod * condition * oneBackFP + 
                                  (1+foreperiod*condition*oneBackFP|ID),
                                data = data2,
                                family = gaussian(),
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.9),
                                warmup = 2000,
                                iter = 12000,
                                cores = -1)


options(digits = 5)
summary(b_one_back_fp)
options(options_defaults)




#===========================================================================================#
#=================================== 9. Accuracy ============================================
#===========================================================================================#

fpaccglmer <- buildmer(acc_result ~ foreperiod * condition + 
                         (1+foreperiod*condition|ID), 
                       data=dataAcc,
                       family = binomial(link = "logit"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite",
                                                         include = ~ foreperiod * condition))


formula(fpaccglmer)

fpaccglmer <- mixed(formula = acc_result ~ foreperiod * condition +
                         (1 | ID),
                       data = dataAcc,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    method = "LRT")


isSingular(fpaccglmer)

summary(fpaccglmer)
anova(fpaccglmer)

