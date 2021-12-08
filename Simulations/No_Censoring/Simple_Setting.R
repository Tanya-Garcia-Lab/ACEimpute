library(impeRfect)
library(lme4)
library(geex)
set.seed(95)

# number of simulations; number of clusters; obs. per cluster
n.sims = 20; n = 1000; m = 3
# vector of fixed effects
beta = c(0, 2, -1)
p = length(beta)
p.X = p - 1
# random error sd; random errors
sigma = 1

data.list = list()
for (i in 1:n.sims) { 
  # when b is NULL, random intercept is generated from N(0, 1)
  b = NULL
  # otherwise, pre-define b as an n-dimensional vector of random intercepts
  # b = rt(n = n, df = 3), for example
  data.list[[i]] = generate_data(n = n, m = m, b = b, p.X = p.X, p.A = 0, beta = beta) 
}

lme.results = ee.results = matrix(data = NA, nrow = n.sims, ncol = p)

for (i in 1:n.sims) {
  print(paste("Sim",i, "starting"))
  sim.data = data.list[[i]]
  # use lmer() to fit linear mixed effects model
  lme.fit = lmer(formula = Y ~ X1 + X2 + (1 | id) - 1, data = sim.data) %>%
    summary()
  lme.results[i, ] = c(coef(lme.fit)[, 1], lme.fit$sigma^2)
  # lme.results[i, 5:7] = sqrt(diag(vcov(lme.fit)))
  
  # use m_estimate() to solve estimating equations defined above
  ee.fit = m_estimate(estFUN = eff_score_vec, 
                      data = sim.data, units = "id",
                      root_control = setup_root_control(start = lme.results[i, ]),
                      outer_args = list(response = "Y", X.names = paste0("X", 1:(p - 1))))
  
  ee.results[i, ] = coef(ee.fit)
  print(paste("Sim",i, "complete"))
  # ee.results[i, 5:7] = sqrt(diag(vcov(ee.fit))[1:3])
}

c(beta, sigma)
colMeans(lme.results)
colMeans(ee.results)

sqrt(vcov(ee.fit))

# mybasis = create_basis(estFUN = eff_score_vec,
#                       data = data.list[[1]], units = "id",
#                       outer_args = list(response = "Y", variant.X = paste0("X", 1:(p - 1))))
# f <- grab_GFUN(create_GFUN(mybasis))
# # Evaluate GFUN at mean and variance: should be close to zero
# f(c(beta, sigma^2))
# f(lme.results[1, ])
# f(ee.results[1, ])
# 
# mybasis = create_basis(estFUN = eff_score_vec,
#                        data = data.list[[2]], units = "id",
#                        outer_args = list(response = "Y", variant.X = paste0("X", 1:(p - 1))))
# f <- grab_GFUN(create_GFUN(mybasis))
# # Evaluate GFUN at mean and variance: should be close to zero
# f(c(beta, sigma^2))
# f(lme.results[2, ])