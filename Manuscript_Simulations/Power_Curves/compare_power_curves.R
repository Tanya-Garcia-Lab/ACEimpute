# clear workspace
rm(list = ls())

setwd("~/Documents/GitHub/ACEimpute/Manuscript_Simulations/Power_Curves")

library(magrittr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

# set seed
set.seed(95)

# sample sizes 5-1200
n <- seq(5, 1200, by = 1)

# estimated slopes on s - T (time to diagnosis) in placebo group
cmi_b1 <- -0.7385014 # CMI-MI
ace_b1 <- -0.7569002 # ACE
cca_b1 <- -1.2580475 # CCA

# assumed variance of delta
varDelta <- 1

# traetment effect = frac*placebo effect
frac = 0.90

# standard normal quantiles
Z_signif <- qnorm(p = 0.05)
Z_power <- qnorm(p = 0.8)

# difference between placebo and treatment groups
cmi_d <- cmi_b1 - frac*cmi_b1
ace_d <- ace_b1 - frac*ace_b1
cca_d <- cca_b1 - frac*cca_b1

# estimated power at sample sizes 5-1200
power_cmib1 <- pnorm(q = (Z_signif - cmi_d * sqrt(n) / sqrt(varDelta)))
power_aceb1 <- pnorm(q = (Z_signif - ace_d * sqrt(n) / sqrt(varDelta)))
power_ccab1 <- pnorm(q = (Z_signif - cca_d * sqrt(n) / sqrt(varDelta)))

# required sample size for power = 80%
req_n_cmi <- n[which.min(abs(power_cmib1 - 0.8))]
req_n_ace <- n[which.min(abs(power_aceb1 - 0.8))]
req_n_cca <- n[which.min(abs(power_ccab1 - 0.8))]

# create x-axis breaks
x_axis_n <- c(seq(0, max(n), by = 200))

# plot power curves
data.frame(Parameter = rep(x = c("CCA", "CMI-MI", "ACE"), each = length(n)),
           SampleSize = rep(n, times = 3),
           Power = c(power_ccab1, power_cmib1, power_aceb1)) %>%
  dplyr::mutate(Method = factor(Parameter, levels = c("CCA", "CMI-MI", "ACE"))) %>%
  ggplot(aes(x = SampleSize, y = Power, linetype = Method, col = Method)) +
  # power curves
  geom_line(size = 1.2) +
  # horizontal line for power = 80%
  geom_hline(yintercept = 0.8, linetype = "solid", color = "gray", size = 0.8) +
  # vertical lines for required sample sizes
  geom_vline(xintercept = req_n_cca, col = "black", linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = req_n_cmi, col = "#E69F00", linetype = "twodash", size = 0.8) +
  geom_vline(xintercept = req_n_ace, col = "#56B4E9", linetype = "solid", size = 0.8) +
  ggthemes::scale_color_colorblind() +
  # set line types
  scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
  # custom y-axis breaks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  # custom x-axis breaks
  scale_x_continuous(limits = c(0, max(n)), breaks = x_axis_n) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  xlab("Sample Size") +
  ylab("Power to Detect\nTreatment Effect") +
  annotate("text", x=800, y=0.30, label = "CCA: n = 391\nCMI-MI: n = 1134\nACE: n = 1079", size = 3.5) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) -> plot_diff_trt_90

# sample sizes 5-600
n <- seq(5, 600, by = 1)

# traetment effect = frac*placebo effect
frac = 0.85

# difference between placebo and treatment groups
cmi_d <- cmi_b1 - frac*cmi_b1
ace_d <- ace_b1 - frac*ace_b1
cca_d <- cca_b1 - frac*cca_b1

# estimated power at sample sizes 5-500
power_cmib1 <- pnorm(q = (Z_signif - cmi_d * sqrt(n) / sqrt(varDelta)))
power_aceb1 <- pnorm(q = (Z_signif - ace_d * sqrt(n) / sqrt(varDelta)))
power_ccab1 <- pnorm(q = (Z_signif - cca_d * sqrt(n) / sqrt(varDelta)))

# required sample size for power = 80%
req_n_cmi <- n[which.min(abs(power_cmib1 - 0.8))]
req_n_ace <- n[which.min(abs(power_aceb1 - 0.8))]
req_n_cca <- n[which.min(abs(power_ccab1 - 0.8))]

# create x-axis breaks
x_axis_n <- c(seq(0, max(n), by = 100))

# plot power curves
data.frame(Parameter = rep(x = c("CCA", "CMI-MI", "ACE"), each = length(n)),
           SampleSize = rep(n, times = 3),
           Power = c(power_ccab1, power_cmib1, power_aceb1)) %>%
  dplyr::mutate(Method = factor(Parameter, levels = c("CCA", "CMI-MI", "ACE"))) %>%
  ggplot(aes(x = SampleSize, y = Power, linetype = Method, col = Method)) +
  # power curves
  geom_line(size = 1.2) +
  # horizontal line for power = 80%
  geom_hline(yintercept = 0.8, linetype = "solid", color = "gray", size = 0.8) +
  # vertical lines for required sample sizes
  geom_vline(xintercept = req_n_cca, col = "black", linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = req_n_cmi, col = "#E69F00", linetype = "twodash", size = 0.8) +
  geom_vline(xintercept = req_n_ace, col = "#56B4E9", linetype = "solid", size = 0.8) +
  ggthemes::scale_color_colorblind() +
  # set line types
  scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
  # custom y-axis breaks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  # custom x-axis breaks
  scale_x_continuous(limits = c(0, max(n)), breaks = x_axis_n) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  xlab("Sample Size") +
  ylab("Power to Detect\nTreatment Effect") +
  annotate("text", x=400, y=0.3, label = "CCA: n = 174\nCMI-MI: n = 504\nACE: n = 480", size = 3.5) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) -> plot_diff_trt_85

# sample sizes 5-300
n <- seq(5, 300, by = 1)

# traetment effect = frac*placebo effect
frac = 0.8

# difference between placebo and treatment groups
cmi_d <- cmi_b1 - frac*cmi_b1
ace_d <- ace_b1 - frac*ace_b1
cca_d <- cca_b1 - frac*cca_b1

# estimated power at sample sizes 5-500
power_cmib1 <- pnorm(q = (Z_signif - cmi_d * sqrt(n) / sqrt(varDelta)))
power_aceb1 <- pnorm(q = (Z_signif - ace_d * sqrt(n) / sqrt(varDelta)))
power_ccab1 <- pnorm(q = (Z_signif - cca_d * sqrt(n) / sqrt(varDelta)))

# required sample size for power = 80%
req_n_cmi <- n[which.min(abs(power_cmib1 - 0.8))]
req_n_ace <- n[which.min(abs(power_aceb1 - 0.8))]
req_n_cca <- n[which.min(abs(power_ccab1 - 0.8))]

# create x-axis breaks
x_axis_n <- c(seq(0, max(n), by = 50))

# plot power curves
data.frame(Parameter = rep(x = c("CCA", "CMI-MI", "ACE"), each = length(n)),
           SampleSize = rep(n, times = 3),
           Power = c(power_ccab1, power_cmib1, power_aceb1)) %>%
  dplyr::mutate(Method = factor(Parameter, levels = c("CCA", "CMI-MI", "ACE"))) %>%
  ggplot(aes(x = SampleSize, y = Power, linetype = Method, col = Method)) +
  # power curves
  geom_line(size = 1.2) +
  # horizontal line for power = 80%
  geom_hline(yintercept = 0.8, linetype = "solid", color = "gray", size = 0.8) +
  # vertical lines for required sample sizes
  geom_vline(xintercept = req_n_cca, col = "black", linetype = "dotted", size = 0.8) +
  geom_vline(xintercept = req_n_cmi, col = "#E69F00", linetype = "twodash", size = 0.8) +
  geom_vline(xintercept = req_n_ace, col = "#56B4E9", linetype = "solid", size = 0.8) +
  ggthemes::scale_color_colorblind() +
  # set line types
  scale_linetype_manual(values=c("dotted", "twodash", "solid")) +
  # custom y-axis breaks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  # custom x-axis breaks
  scale_x_continuous(limits = c(0, max(n)), breaks = x_axis_n) +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  xlab("Sample Size") +
  ylab("Power to Detect\nTreatment Effect") +
  annotate("text", x=200, y=0.30, label = "CCA: n = 98\nCMI-MI: n = 283\nACE: n = 270", size = 3.5) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) -> plot_diff_trt_80

# arrange figures in column
figure = ggarrange(plot_diff_trt_80, 
                   plot_diff_trt_85, 
                   plot_diff_trt_90, 
                   # use common legend
                   common.legend = TRUE, legend = "top",
                   # add labels to indicate treatment effect, offset to avoid overlapping other figure text
                   labels = c("20%", "15%", "10%"), label.y = 1.1, label.x = 0.9, ncol = 1)

# add x- and y-axis labels
figure = annotate_figure(figure, 
                         left = text_grob("Power to detect treatment effect", rot = 90),
                         bottom = text_grob("Required sample size"))

# save annotated figure
ggsave(filename = "compare_power_curves.png", plot = figure, device = "png")