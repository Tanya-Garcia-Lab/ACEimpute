library(tidyverse)
library(latex2exp)

setwd("~/Documents/GitHub/Random-Error-Imputation/Simulations/Censoring/Impute_T/Results")

heavy_cens = read.csv("IT_PhRMASetting_HeavyCens_31Jan22.csv")

light_cens = read.csv("IT_PhRMASetting_LigthCens_18Feb22.csv")

heavy_cens_long = heavy_cens %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

heavy_cens_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

light_cens_long = light_cens %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Light Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

light_cens_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

long_data = rbind(heavy_cens_long, light_cens_long)

# Replicate Figure S2 from the Supplemental Materials for Lotspeich, Grosser and Garcia (2022+)
alpha_plot = long_data %>%
  filter(param == "alpha") %>%
  ggplot(aes(y = est, x = setting, fill = method)) +
  geom_boxplot() + geom_hline(yintercept = 1, linetype = 2, color = "black") +
  theme_bw() +
  scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
  theme(legend.position = "top", axis.title.x = element_blank()) +
  ylab(TeX("Parameter estimate: $\\hat{\\alpha}$"))

alpha_plot

beta_plot = long_data %>%
  filter(param == "beta1") %>%
  ggplot(aes(y = est, x = setting, fill = method)) +
  geom_boxplot() + geom_hline(yintercept = 0.25, linetype = 2, color = "black") +
  theme_bw() +
  scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
  theme(legend.position = "top", axis.title.x = element_blank()) +
  ylab(TeX("Parameter estimate: $\\hat{\\beta}$"))

beta_plot

sigma2_plot = long_data %>%
  filter(param == "sigma2") %>%
  ggplot(aes(y = est, x = setting, fill = method)) +
  geom_boxplot() + geom_hline(yintercept = 1, linetype = 2, color = "black") +
  theme_bw() +
  scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
  theme(legend.position = "top", axis.title.x = element_blank()) +
  ylab(TeX("Parameter estimate: $\\hat{\\sigma}^2$"))

sigma2_plot