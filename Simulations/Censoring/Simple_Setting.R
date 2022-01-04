# need this for simulation of Cox model outcome
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
library(lme4)
library(geex)


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

# create time_since_event variable, age - t
long.data$time_since_event = with(long.data, age - t)

# inspect data frame
head(long.data)

# check censoring rate
1 - mean(long.data$delta)

# sanity check with lm()
# truth = (0.5, 2, -1, 1.5)
# NOTE: I use the covarirate (age - t), which correctly follows the DGM
lm(formula = y ~ z_y1 + z_y2 + time_since_event, data = long.data) %>% coef()

# sanity check with coxph()
# truth = (1, -0.5)
survival::coxph(formula = survival::Surv(t, delta) ~ z_t1 + z_t2, data = long.data) %>% coef()

## HI SARAH! :)
## NEXT STEPS:
# impute using imputeCensoRd::condl_mean_impute()
# analyze with impeRfect::eff_score_vector()


# get initial parameter estimates with lmer()
lme.fit = lmer(formula = y ~ z_y1 + z_y2 + time_since_event + (1 | id) - 1, data = long.data)
lme.sum = summary(lme.fit)
lme.est = c(coef(lme.sum)[, 1], lme.sum$sigma^2)
lme.est

# use m_estimate() to solve estimating equations defined above
ee.fit = m_estimate(estFUN = eff_score_vec, 
                    data = long.data, units = "id",
                    root_control = setup_root_control(start = lme.est),
                    outer_args = list(response = "y", 
                                      X.names = c("z_y1", "z_y2", "time_since_event"),
                                      cens = "delta"))

# truth = c(2, -1, 1.5, 1)
coef(ee.fit)

