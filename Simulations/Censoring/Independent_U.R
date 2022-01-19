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

# imputation error U
# is this SD too big?
U = rnorm(n = n, mean = 0, sd = 1)
long.data$U = (1 - long.data$delta)*rep(U, each = m)

# if censored, add measurement error
long.data$t_star = with(long.data, t + U)
long.data$time_to_event_star = with(long.data, time_to_event - U)

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

# get initial parameter estimates with lmer()
lme.fit <- lmer(formula = y ~ z_y1 + z_y2 + time_to_event_star + (1 | id) - 1, data = long.data)
lme.sum <- summary(lme.fit)
lme.est <- c(coef(lme.sum)[, 1], lme.sum$sigma^2)
lme.est

# use m_estimate() to solve estimating equations defined above
ee.fit <- m_estimate(estFUN = eff_score_vec,
                     data = long.data, units = "id",
                     root_control = setup_root_control(start = lme.est),
                     outer_args = list(response = "y",
                                       X.names = c("z_y1", "z_y2", "time_to_event_star"),
                                       # This line is the only real change needed for censored data
                                       cens = "delta"))

# truth = c(2, -1, 1, 1)
coef(ee.fit)

# empirical sandwich estimator of the standard errors
sqrt(diag(vcov(ee.fit)))