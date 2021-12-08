
# Run once: install.packages("devtools")
devtools::install_github(repo = "Tanya-Garcia-Lab/Imputing-Censored-Covariates", subdir = "imputeCensoRd")
library(tidyverse)

set.seed(114)
# number of subjects
n = 1000
# nimber of visits per subject
m = 3

## STEP 1: generate z.t, time-independent covariates
# log hazard ratio for simulation of t, the survival outcome
logHR <- c(1, -0.5)
# number of covariates in z.t
p.logHR <- length(logHR)
# simulate z.t as iid from Normal(0, 1)
z.t <- rnorm(n = n*p.logHR, mean = 0, sd = 1) %>%
  matrix(ncol = p.logHR)
colnames(z.t) = c("z_t1", "z_t2")

## STEP 2: generate t from Cox simulation
t <- imputeCensoRd::cox_simulation(n = n, logHR = logHR, A = z.t, dist = "Exponential", lambda = 5)

# Sanity check on Cox simulation
# survival::coxph(formula = survival::Surv(t, rep(1, n)) ~ z.t)

## STEP 3: generate y, the longitudinal response, and z.y, the time-dependent covariates
# intercept and coefficients on z.y
beta = c(0.5, 2, -1)
p.beta = length(beta)
long.data = impeRfect::generate_data(n = n, m = m, p.X = p.beta - 1, p.A = 0, beta = beta)

# sanity check on linear model simulation
# lm(formula = Y ~ X1 + X2, data = long.data)

# coefficient on t
beta.t = 1.5
# add beta.t * t to Y
long.data = long.data %>%
  mutate(t = rep(t, each = m)) %>%
  mutate(Y = Y + beta.t*t)

# add time-independent covariates to data
long.data = long.data %>%
  cbind(kronecker(z.t, rep(1, 3)))

# renaming columns to match Overleaf documentation
colnames(long.data) = c("id", "y", paste0("z_y", 1:(p.beta - 1)),
                        "t", paste0("z_t", 1:p.logHR))

## STEP 3: generate c, the censoring mechanism
c <- runif(n = n, min = 0, max = 1)
delta <- as.numeric(t <= c)
w <- ifelse(delta, t, c)

# add censoring information to data
long.data = long.data %>%
  mutate(c = rep(c, each = m),
         delta = rep(delta, each = m),
         w = rep(w, each = m))

# inspect data frame
head(long.data)
# censoring rate
1 - mean(long.data$delta)
# another sanity check
# truth = (0.5, 2, -1, 1.5)
lm(formula = y ~ z_y1 + z_y2 + t, data = long.data) %>% coef()