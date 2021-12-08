#' Calculate the efficient score function for theta
#'
#' Take in in dataframe 'data', which corresponds to data from a
#' single cluster, then return the efficient score function for
#' theta = (beta, alpha, sigma^2) corresponding to that cluster.
#' Intended for use with m_estimate() function
#' Right now, intended for use in linear mixed model when we do
#' not want to specify a distribution for the random effects.
#' Also, at the moment, this function only allows for random
#' intercept.
#'
#' @param data dataframe for one cluster of observations
#' @param response character denoting the column name for the model response
#' @param X.names string of characters denoting the columns for the model covariates that are associated with fixed effects
#' @param cens character denoting the column name for the censoring indicator; 1 is not censored and 0 is censored
# @param Z.names string of characters denoting the columns for the model covariates that are associated with random effects
#'
#' @importFrom MASS ginv
#'
#' @export
eff_score_vec = function(data, response, X.names, cens = NULL) {
                         # Z.names = NULL) {
  # number of observations in cluster
  m = nrow(data)

  # cens.bool = 1 if cens != NULL and data[1, cens] == 0
  cens.bool = !is.null(cens)
  if (cens.bool) { cens.bool = (data[1, cens] == 0) }

  # number of fixed effects
  p = length(X.names)
  # number of random effects
  q = 1 #+ length(Z.names)

  y = as.matrix(data[, response])
  # fixed effects design matrix
  X = as.matrix(data[, X.names])

  # random effects design matrix
  Z = matrix(1, m, 1)
  if (cens.bool) {
    Z = cbind(1, Z)
    q = q + 1
  }
  # if(!is.null(Z.names)) { Z = cbind(Z, as.matrix(data[, Z.names])) }

  ## THIS IS HOW I HAD TO WRITE THE FUNCTION FOR IT BE USED WITH
  ## GEEX - THE FUNCTION ARGUMENT IS THETA AND IT RETURNS A VECTOR
  ## THAT IS THE SUMMAND FOR THE ESTIMATING EQUATION OF INTEREST
  ## IT'S KINDA WEIRD BUT I GUESS IT MAKES SENSE
  function(theta) {
    # variance of random error
    sigma2 = theta[p+1]

    # transform Y; explained in Overleaf
    # w = t(Z) %*% y / sigma2

    # calculate zeta, the fixed effects component of linear predictor
    zeta = X %*% theta[1:p]

    # conditional expectation of y given w, x, s, t, z
    cond.expect.y <- zeta + Z %*% MASS::ginv(t(Z) %*% Z) %*% t(Z) %*% (y -  zeta)

    # conditional expectation of y'y given w, x, s, t, z
    cond.expect.y2 <- sigma2*(m - q) + sum(cond.expect.y^2)

    ## BETA
    # efficient score vector for beta
    eff.score.beta <- t(X) %*% (y - cond.expect.y)/sigma2

    ## SIGMA
    # efficient score for sigma^2
    eff.score.sigma2 <- (0.5*(sum(y^2) - cond.expect.y2) - t(zeta) %*% (y - cond.expect.y))/(sigma2^2)

    c(eff.score.beta, eff.score.sigma2)
  }
}
