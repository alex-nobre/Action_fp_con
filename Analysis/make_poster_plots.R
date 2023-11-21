
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
library(ggsignif)
library(viridis)

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

# 

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


RT_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 4.3, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 4.1, width = 0.1, geom = "errorbar") + 
  labs(title = "Experiment 3 (n = 28)",
       x = "",
       y = "",
       color = "Condition") +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

# Save
ggplot2::ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/RT_by_condition_exp3.tiff",
                RT_by_condition,
                width = 25,
                height = 16.66,
                units = "cm")



#=============================================================================================#
#====================================== Plots em português ===================================#
#=============================================================================================#

RT_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.6, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.4, width = 0.1, geom = "errorbar") + 
  labs(title = "TR",
       x = "",
       y = "TR médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))

# Save
ggplot2::ggsave("G:/My Drive/Post-doc/Orientacoes/Sabricia/PDPD_2022-2023/Simposio_pos/RT_por_condicao.tiff",
                RT_by_condition,
                width = 25,
                height = 16.66,
                units = "cm")

error_by_condition <- ggplot(data = summaryDataAcc %>%
                               group_by(ID, foreperiod, condition) %>%
                               summarise(errorRate = mean(errorRate)),
                             aes(x = foreperiod,
                                 y = errorRate,
                                 color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.6, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.4, width = 0.1, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Taxa de erro média",
       color = "Condição",
       title = "Erros") +
  scale_color_manual(values = c("orange", "blue"),
                     label = c("Externa", "Ação")) +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "bottom")

cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(error_by_condition + 
                                                          theme(legend.title = element_text(size = rel(1.6)),
                                                                legend.text = element_text(size = rel(1.4))))), "guide-box")

xaxis_title <- text_grob(error_by_condition$labels$x,
                         just = "top",
                         size = (error_by_condition + 
                                   theme(axis.title = element_text(size = rel(1.4))))$theme$axis.title$size * 11) # 11 is the base size in theme_grey


xaxis_title_margin <- unit(2, "line")

rt_error_plots <- arrangeGrob(arrangeGrob(RT_by_condition + theme(legend.position = "none",
                                                                  axis.text = element_text(size = rel(1.2)),
                                                                  axis.title = element_text(size = rel(1.4)),
                                                                  plot.title = element_text(size = rel(1.5)),
                                                                  axis.title.x = element_blank()),
                                          error_by_condition + theme(legend.position = "none",
                                                                     axis.text = element_text(size = rel(1.2)),
                                                                     axis.title = element_text(size = rel(1.4)),
                                                                     plot.title = element_text(size = rel(1.5)),
                                                                     axis.title.x = element_blank()) + labs(x = ""),
                                          nrow = 1,
                                          widths = c(4/8, 4/8)),
                              xaxis_title,
                              cond_legend,
                              heights = unit.c(unit(1, "null"),
                                               grobHeight(xaxis_title) + xaxis_title_margin, grobHeight(xaxis_title) + xaxis_title_margin),
                              nrow = 3)


# Save
ggplot2::ggsave("G:/My Drive/Post-doc/Orientacoes/Sabricia/PDPD_2022-2023/Simposio_pos/RT_erros.tiff",
                rt_error_plots,
                width = 25,
                height = 16.66,
                units = "cm")

# Check for influence of external fixation duration
extDurCheck <- ggplot(data=filter(data,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(1.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.6)),
        legend.text = element_text(size = rel(1.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "none") +
  labs(x = "Duração da fixação na condição externa (s)",
       y = "TR (s)",
       color = "FP") +
  facet_wrap(~foreperiod,
             labeller = as_labeller(c(`1` = "FP = 1.0 s",
                                      `2.8` = "FP = 2.8 s"))) +
  scale_color_viridis(discrete = TRUE, begin = 0.20, end = 0.85)
ggsave("G:/My Drive/Post-doc/Orientacoes/Sabricia/PDPD_2022-2023/Simposio_pos/extfixduration.tiff",
       extDurCheck,
       width = 25,
       height = 16.66,
       units = "cm")

# Check for influence of latency of action key press on RT
actionTrigPress <- ggplot(data=filter(data,condition=='action',action_trigger.rt < 5.0),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(1.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.3)),
        legend.title = element_text(size = rel(1.6)),
        legend.text = element_text(size = rel(1.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "none") +
  labs(x = "Tempo para início voluntário da tentativa (s)",
       y = "TR (s)",
       color = "FP") +
  facet_wrap(~foreperiod,
             labeller = as_labeller(c(`1` = "FP = 1.0 s",
                                      `2.8` = "FP = 2.8 s"))) +
  scale_color_viridis(discrete = TRUE, begin = 0.20, end = 0.85)
  
ggsave("G:/My Drive/Post-doc/Orientacoes/Sabricia/PDPD_2022-2023/Simposio_IC/action_trig_press.tiff",
       actionTrigPress,
       width = 25,
       height = 16.66,
       units = "cm")
