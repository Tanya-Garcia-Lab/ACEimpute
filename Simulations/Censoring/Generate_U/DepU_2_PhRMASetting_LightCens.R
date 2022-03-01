libFolder = "/nas/longleaf/home/lotspeic/R/"
library(impeRfect, lib.loc = libFolder)
library(lme4)
library(geex, lib.loc = libFolder)
library(magrittr)

args <- commandArgs(TRUE)
sim_seed <- 918 + as.integer(args)
set.seed(sim_seed)

# number of subjects
n <- 1000
# number of visits per subject
m <- 3

## REGRESSION PARAMETERS
# log hazard ratio on z.t for simulation of t, the survival outcome
logHR <- c(1, -0.5)
rate.t <- 0.5
rate.c <- 1 # rate.c = 1 --> ~20% censoring, rate.c = 5 --> ~45% censoring
# intercept and coefficients on z.y (in that order)
beta <- c(1, 0.25)
# coefficient on t for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error
sigma <- 1

num_sim <- 1000
save_res <- data.frame(sim = 1:num_sim, censored = NA,
                       beta1_ee = NA, se_beta1_ee = NA, alpha_ee = NA, se_alpha_ee = NA, sigma2_ee = NA,
                       beta1_lme = NA, se_beta1_lme = NA, alpha_lme = NA, se_alpha_lme = NA, sigma2_lme = NA)

data_list = list()
for (s in 1:num_sim) {
  ## GENERATE LONGITUDINDAL DATA USING impeRfect::DGM_2
  data_list[[s]] <- DGM_2(n = n, m = m, b = NULL,
                          beta = beta, alpha = alpha, sigma = sigma,
                          logHR = logHR, rate.t = rate.t, rate.c = rate.c)
}

for (s in 1:num_sim) {
  long.data = data_list[[s]]
  
  # imputation error U
  # is this SD too big?
  unique_t = long.data[which(long.data$visit == 1), "t"]
  U = rnorm(n = n, mean = 0, sd = unique_t)
  long.data$U = (1 - long.data$delta)*rep(U, each = m)
  
  # if censored, add measurement error
  long.data$imp = with(long.data, t + U)
  long.data$time_to_event_imp = with(long.data, (visit - 1) - imp)
  
  # check censoring rate
  save_res[s, "censored"] <- 1 - mean(long.data$delta)
  
  # analyze with impeRfect::eff_score_vector(), as done in the following code
  # get initial parameter estimates with lmer()
  lme.fit <- lmer(formula = y ~ z_y1 + time_to_event_imp + (1 | id) - 1, data = long.data)
  lme.sum <- summary(lme.fit)
  lme.est <- save_res[s, c("beta1_lme", "alpha_lme", "sigma2_lme")] <- c(coef(lme.sum)[, 1], lme.sum$sigma^2)
  save_res[s, c("se_beta1_lme", "se_alpha_lme")] <- sqrt(diag(vcov(lme.fit)))
  
  # use m_estimate() to solve estimating equations defined above
  ee.fit <- m_estimate(estFUN = eff_score_vec,
                       data = long.data, units = "id",
                       root_control = setup_root_control(start = lme.est),
                       outer_args = list(response = "y",
                                         X.names = c("z_y1", "time_to_event_imp"),
                                         # This line is the only real change needed for censored data
                                         cens = "delta"))
  save_res[s, c("beta1_ee", "alpha_ee", "sigma2_ee")] <- coef(ee.fit)
  save_res[s, c("se_beta1_ee", "se_alpha_ee")] <- sqrt(diag(vcov(ee.fit)[-3, -3]))
  
  if (s %% 25 == 0) print(paste("Simulation", s, "complete!"))
  
  write.csv(save_res, "~kylegrosser/Documents/GitHub/Random-Error-Imputation/Simulations/Censoring/Generate_U/Results/DepU_2_PhRMASetting_LightCens.csv", row.names = F)
}

library(magrittr)
save_res %>%
  dplyr::select(-sim) %>%
  dplyr::summarize_all(mean)