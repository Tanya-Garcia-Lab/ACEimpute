# need this for simulation of Cox model outcome
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
library(tidyverse)

set.seed(114)
# number of subjects
n = 1000
# number of visits per subject
m = 3

## REGRESSION PARAMETERS
# log hazard ratio on z.t for simulation of t, the survival outcome
logHR <- c(1, -0.5)
# intercept and coefficients on z.y (in that order)
beta = c(0.5, 2, -1)
# coefficient on t for simulation of y
alpha = 1.5
# standard deviation of epsilon, the random error 
sigma = 1

## GENERATE LONGITUDINDAL DATA USING impeRfect::DGM_2
long.data = DGM_2(n = n, m = m, b = NULL, 
                  logHR = logHR, beta = beta, 
                  alpha = alpha, sigma = sigma)

# inspect data frame
head(long.data)

# check censoring rate
1 - mean(long.data$delta)

# sanity check with lm()
# truth = (0.5, 2, -1, 1.5)
# NOTE: I use the covarirate (age - t), which correctly follows the DGM
lm(formula = y ~ z_y1 + z_y2 + I(age - t), data = long.data) %>% coef()

# sanity check with coxph()
# truth = (1, -0.5)
survival::coxph(formula = survival::Surv(t, delta) ~ z_t1 + z_t2, data = long.data) %>% coef()

## HI SARAH! :)
## NEXT STEPS:
# impute using imputeCensoRd::condl_mean_impute()
# analyze with impeRfect::eff_score_vector()