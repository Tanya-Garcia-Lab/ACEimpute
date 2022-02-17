# need this for simulation of Cox model outcome
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
library(lme4)
library(geex)

set.seed(114)
# number of subjects
n <- 1000
# number of visits per subject
m <- 3

## REGRESSION PARAMETERS
# log hazard ratio on z.t for simulation of t, the survival outcome
logHR <- c(1, -0.5)
rate.t <- 0.5
rate.c <- 1
# intercept and coefficients on z.y (in that order)
beta <- c(0, 2, -1)
# coefficient on t for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error 
sigma <- 1

## GENERATE LONGITUDINDAL DATA USING impeRfect::DGM_2
long.data <- DGM_2(n = n, m = m, b = NULL, 
                   beta = beta, alpha = alpha, sigma = sigma,
                   logHR = logHR, rate.t = rate.t, rate.c = rate.c)

# inspect data frame
head(long.data)

# what does the event time t look like? is it realistic?
summary(long.data$t)

# check censoring rate
1 - mean(long.data$delta)

# sanity check with lm()
# truth = c(0, 2, -1, 1)
summary(lmer(formula = y ~ z_y1 + z_y2 + time_to_event + (1 | id), data = long.data))$coefficients[, 1]
# truth = 1
summary(lmer(formula = y ~ z_y1 + z_y2 + time_to_event + (1 | id), data = long.data))$sigma

# sanity check with coxph()
# truth = c(1, -0.5)
coef(survival::coxph(formula = survival::Surv(w, delta) ~ z_t1 + z_t2, data = long.data))

# impute using imputeCensoRd::condl_mean_impute()
## fit imputation model for W|Z_t1, Z_t2
imp_mod <- survival::coxph(formula = survival::Surv(w, delta) ~ z_t1 + z_t2, data = long.data)
## Impute censored covariate t
imp_data <- imputeCensoRd::condl_mean_impute(fit = imp_mod, obs = "w", event = "delta", addl_covar = c("z_t1", "z_t2"), data = long.data)
## recalculate time_to_event using imputed value ("imp")
imp_data$time_to_event_imp <- with(imp_data, (visit - 1) - imp)

# analyze with impeRfect::eff_score_vector()
## get initial parameter estimates using imputed time_to_event_imp with lmer()
lme.fit <- lmer(formula = y ~ z_y1 + z_y2 + time_to_event_imp + (1 | id) - 1, data = imp_data)
lme.sum <- summary(lme.fit)
lme.est <- c(coef(lme.sum)[, 1], lme.sum$sigma^2)
# truth = c(2, -1, 1, 1)
lme.est

# use m_estimate() to solve estimating equations defined above
ee.fit <- m_estimate(estFUN = eff_score_vec,
                     data = imp_data, units = "id",
                     root_control = setup_root_control(start = lme.est),
                     outer_args = list(response = "y",
                                       X.names = c("z_y1", "z_y2", "time_to_event_imp"),
                                       # This line is the only real change needed for censored data
                                       cens = "delta"))

# truth = c(2, -1, 1, 1)
coef(ee.fit)

# empirical sandwich estimator of the standard errors
sqrt(diag(vcov(ee.fit)))