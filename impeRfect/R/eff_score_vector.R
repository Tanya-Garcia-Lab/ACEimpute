#' This function calculates the efficient score vector proposed by
#' Grosser et al. (2022)
#'
#' Estimates the parameters of a linear mixed model when we do
#' not want to specify a distribution for the random effects and
#' when one covariate is (possibly) subject to random right-censoring
#' Takes in dataframe 'data', which corresponds to data from a
#' single cluster, then returns the value of efficient score vector for
#' theta = (beta, sigma^2) corresponding to that cluster.
#' Intended for use with m_estimate() function
#' At the moment, this function only allows for random intercept.
#'
#' @param data dataframe for one cluster of observations
#' @param response character denoting the column name for the model response
#' @param X.names string of characters denoting the columns for the model covariates that are associated with fixed effects
#' @param cens character denoting the column name for the censoring indicator, which is 1 if not censored and 0 if censored.
#' Default value is NULL, in which case the cluster is considered not censored.
#' @param Z.names string of characters denoting the columns for the model covariates that are associated with random effects.
#' Default value is NULL, in which case the random effects design matrix is assumed to be a column of 1's
#'
#' @importFrom MASS ginv
#'
#' @export
eff_score_vec = function(data, response,
                         X.names, Z.names = NULL, cens = NULL) {
  # number of observations in cluster
  m = nrow(data)

  # cens.bool = 1 if cens != NULL and data[1, cens] == 0
  cens.bool = !is.null(cens)
  if (cens.bool) { cens.bool = (data[1, cens] == 0) }

  # number of fixed effects
  p = length(X.names)
  # number of random effects
  q = 1 + length(Z.names)

  # longitudinal response
  y = as.matrix(data[, response])
  # fixed effects design matrix
  X = as.matrix(data[, X.names])

  # random effects design matrix
  if (is.null(Z.names)) {
    # if Z.names not provided, Z = 1_m
    Z = matrix(1, m, 1)
  }
  else {
    # otherwise, use the design matrix specified by Z.names
    Z = as.matrix(data[, Z.names])
  }
  # add leading column of 1_m if censored
  if (cens.bool) {
    Z = cbind(1, Z)
    q = q + 1
  }

  # The following lines are written so that this function is
  #  compatible with the function m_estimate() from the geex package
  function(theta) {
    # variance of random error
    sigma2 = theta[p+1]

    # calculate zeta, the fixed effects component of linear predictor
    zeta = X %*% theta[1:p]

    # conditional expectation of y given w, x, s, t, z
    cond.expect.y <- zeta + Z %*% MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% (y -  zeta)

    # conditional expectation of y'y given w, x, s, t, z
    cond.expect.y2 <- sigma2*(m - q) + sum(cond.expect.y^2)

    # efficient score vector for beta
    eff.score.beta <- t(X) %*% (y - cond.expect.y)/sigma2

    # efficient score for sigma^2
    eff.score.sigma2 <- (0.5*(sum(y^2) - cond.expect.y2) - t(zeta) %*% (y - cond.expect.y))/(sigma2^2)

    # return efficient vector for theta = c(beta, sigma^2)
    c(eff.score.beta, eff.score.sigma2)
  }
}
