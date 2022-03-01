library(tidyverse)
library(latex2exp)

setwd("~/Documents/GitHub/Random-Error-Imputation/Simulations/Censoring/Generate_U/Results")

Dep1_heavy = read.csv("DepU_1_PhRMASetting_HeavyCens_21Feb22.csv")
Dep1_light = read.csv("DepU_1_PhRMASetting_LightCens_21Feb22.csv")
Dep2_heavy = read.csv("DepU_2_PhRMASetting_HeavyCens_21Feb22.csv")
IndU_heavy = read.csv("IndU_PhRMASetting_HeavyCens_21Feb22.csv")
IndU_light = read.csv("IndU_PhRMASetting_LightCens_21Feb22.csv")

Dep1_heavy_long = Dep1_heavy %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Depend. 1, Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

Dep1_heavy_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

Dep1_light_long = Dep1_light %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Depend. 1, Light Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

Dep1_light_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

Dep2_heavy_long = Dep2_heavy %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Depend. 2, Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

Dep2_heavy_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

IndU_heavy_long = IndU_heavy %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Indep. U, Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

IndU_heavy_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

IndU_heavy_long = IndU_heavy %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Indep. U, Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

IndU_heavy_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

IndU_light_long = IndU_light %>%
  gather(key = "param0", value = "est", c("beta1_ee", "alpha_ee", "sigma2_ee", 
                                          "beta1_lme", "alpha_lme", "sigma2_lme")) %>%
  mutate(method = rep(c("ee", "reml"), each = 3000),
         param = rep(c("beta1", "alpha", "sigma2"), each = 1000, 2),
         setting = "Indep. U, Heavy Censoring") %>%
  select(sim, setting, censored, method, param, est) %>%
  mutate(method = factor(method, levels = c("ee", "reml"),
                         labels = c("Random Error\nImputation", 
                                    "Conditional Mean\nImputation")),
         param = factor(param, levels = c("beta1", "alpha", "sigma2")))

IndU_light_long %>%
  group_by(method, param) %>%
  summarize(mean = mean(est))

# alpha_plot = long_data %>%
#   filter(param == "alpha") %>%
#   ggplot(aes(y = est, x = setting, fill = method)) +
#   geom_boxplot() + geom_hline(yintercept = 1, linetype = 2, color = "black") +
#   theme_bw() +
#   scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
#   theme(legend.position = "top", axis.title.x = element_blank()) +
#   ylab(TeX("Parameter estimate: $\\hat{\\alpha}$"))
# 
# alpha_plot
# 
# beta_plot = long_data %>%
#   filter(param == "beta1") %>%
#   ggplot(aes(y = est, x = setting, fill = method)) +
#   geom_boxplot() + geom_hline(yintercept = 0.25, linetype = 2, color = "black") +
#   theme_bw() +
#   scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
#   theme(legend.position = "top", axis.title.x = element_blank()) +
#   ylab(TeX("Parameter estimate: $\\hat{\\beta}$"))
# 
# beta_plot
# 
# sigma2_plot = long_data %>%
#   filter(param == "sigma2") %>%
#   ggplot(aes(y = est, x = setting, fill = method)) +
#   geom_boxplot() + geom_hline(yintercept = 1, linetype = 2, color = "black") +
#   theme_bw() +
#   scale_fill_manual(values = rcartocolor::carto_pal(n = 4, name = "Safe"), name = "Method") +
#   theme(legend.position = "top", axis.title.x = element_blank()) +
#   ylab(TeX("Parameter estimate: $\\hat{\\sigma}^2$"))
# 
# sigma2_plot