# need this for simulation of Cox model outcome
# devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
# library(lme4)
# library(geex)

set.seed(114)
n = 400

## REGRESSION PARAMETERS
# log hazard ratio on z.t for simulation of t, the survival outcome
logHR <- c(1, -0.5)
rate.t <- 0.5
# rate.c = 1 gives ~20% censoring
# rate.c = 3 gives ~40% censoring
# rate.c = 10 gives ~60% censoring
rate.c <- 10
# intercept and coefficients on z.y (in that order)
beta <- c(0, 2, -1)
# coefficient on t for simulation of y
alpha <- 1
# standard deviation of epsilon, the random error 
sigma <- 1

## GENERATE LONGITUDINDAL DATA USING impeRfect::DGM_2
long.data.m3 <- DGM_2(n = 133, m = 3, b = NULL, 
                      beta = beta, alpha = alpha, sigma = sigma,
                      logHR = logHR, rate.t = rate.t, rate.c = rate.c)
long.data.m3$id = rep(1:133, each = 3)

long.data.m4 <- DGM_2(n = 133, m = 4, b = NULL, 
                      beta = beta, alpha = alpha, sigma = sigma,
                      logHR = logHR, rate.t = rate.t, rate.c = rate.c)
long.data.m4$id = rep(134:266, each = 4)

long.data.m5 <- DGM_2(n = 134, m = 5, b = NULL, 
                      beta = beta, alpha = alpha, sigma = sigma,
                      logHR = logHR, rate.t = rate.t, rate.c = rate.c)
long.data.m5$id = rep(267:400, each = 5)

long.data = rbind(long.data.m3, long.data.m4, long.data.m5)

# inspect data frame
head(long.data)

# what does the event time t look like? is it realistic?
summary(long.data$t)

# check censoring rate
1 - mean(long.data$delta)

# sanity check with lmer()
# truth = c(2, -1, 1)
summary(lmer(formula = y ~ z_y1 + z_y2 + time_to_event + (1 | id) - 1, data = long.data))$coefficients[, 1]
# truth = 1
summary(lmer(formula = y ~ z_y1 + z_y2 + time_to_event + (1 | id) - 1, data = long.data))$sigma

# sanity check with coxph()
# truth = c(1, -0.5)
invariant.data <- long.data[which(long.data$visit == 1), ]
coef(survival::coxph(formula = survival::Surv(w, delta) ~ z_t1 + z_t2, data = invariant.data))

invariant.data$U = rnorm(n = n, mean = 0, sd = invariant.data$t)
invariant.data$U = with(invariant.data, (1 - delta)*U)
invariant.data = invariant.data[, c("id", "U")]

new_data <- dplyr::left_join(long.data, invariant.data, by = c("id"))
new_data$imp <- with(new_data, t + U)
new_data$time_to_event_imp <- with(new_data, (visit - 1) - imp)

# analyze with impeRfect::eff_score_vector()
## get initial parameter estimates using imputed time_to_event_imp with lmer()
lme.fit <- lmer(formula = y ~ z_y1 + z_y2 + time_to_event_imp + (1 | id) - 1, data = new_data)
lme.sum <- summary(lme.fit)
lme.est <- c(coef(lme.sum)[, 1], lme.sum$sigma^2)
# truth = c(2, -1, 1, 1)
lme.est

# use m_estimate() to solve estimating equations defined above
ee.fit <- m_estimate(estFUN = eff_score_vec,
                     data = new_data, units = "id",
                     root_control = setup_root_control(start = lme.est),
                     outer_args = list(response = "y",
                                       X.names = c("z_y1", "z_y2", "time_to_event_imp"),
                                       # This line is the only real change needed for censored data
                                       cens = "delta"))

# truth = c(2, -1, 1, 1)
coef(ee.fit)

# empirical sandwich estimator of the standard errors
sqrt(diag(vcov(ee.fit)))