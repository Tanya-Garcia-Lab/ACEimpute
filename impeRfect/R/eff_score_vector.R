#' Calculate the efficient score function for theta
#'
#' Take in in dataframe 'data', which corresponds to data from a
#' single cluster, then return the efficient score function thot
#' beta corresponding to that cluster.
#' Intended for use with geex() function
#'
#' @param data dataframe for one cluster of observations
#' @param response character denoting the column name for the model response
#' @param invariant.X character denoting the column name for the time-invariant covariate
#' @param variant.X character denoting the column name for the time-variant covariate
# @param X.names string of characters denoting the columns for the model covariates that are associated with fixed effects
# @param Z.names string of characters denoting the columns for the model covariates that are associated with random effects
#'
#' @export
eff_score_vec = function(data, response, invariant.X = NULL, variant.X) {
                         #, X.names, Z.names = NULL) {
  # number of observations in cluster
  m = nrow(data)

  # is invariant.X null?
  invariant.bool = is.null(invariant.X)

  # number of fixed effects
  p = 1 + length(invariant.X) + length(variant.X)
  # number of random effects
  q = 1 #+ length(Z.names)

  y = as.matrix(data[, response])
  # fixed effects design matrix
  X = as.matrix(data[, variant.X])

  # random effects design matrix
  Z = matrix(1, m, 1)
  # if(!is.null(Z.names)) { Z = cbind(Z, as.matrix(data[, Z.names])) }
  # transpose of Z
  Z = t(Z)
  # columns (q + 1) through m of t(Z)
  Zr = matrix(Z[, (q + 1):m], q)

  ## THIS IS HOW I HAD TO WRITE THE FUNCTION FOR IT BE USED WITH
  ## GEEX - THE FUNCTION ARGUMENT IS THETA AND IT RETURNS A VECTOR
  ## THAT IS THE SUMMAND FOR THE ESTIMATING EQUATION OF INTEREST
  ## IT'S KINDA WEIRD BUT I GUESS IT MAKES SENSE
  function(theta) {
    # variance of random error
    sigma2 = theta[p+1]

    # transform Y
    w = Z %*% y / sigma2
    v = y[(q+1):m] / sigma2

    ## BETA
    # calculate zeta, the fixed effects component of linear predictor
    zeta = theta[1] + X %*% theta[(2 + !invariant.bool):p]
    if (!invariant.bool) { zeta = zeta + data[1, invariant.X] * theta[2] }

    # conditional expectation of y given w, x, s, t, z, b
    cond.expect.y <- zeta + t(Z) %*% solve(Z %*% t(Z)) %*% (sigma2 * w - Z %*% zeta)
    # conditional expectation of y'y given w, x, s, t, z, b
    cond.expect.sum.y2 <- sigma2*(m - q) + sum(cond.expect.y^2)

    ## BETA
    # intercept
    eff.score.beta0 <- sum(y - cond.expect.y)
    # fixed effect for time-invariant covariate
    if (invariant.bool) { eff.score.beta1 = NULL }
    else { eff.score.beta1 <- data[1, invariant.X] * eff.score.beta0 }
    # fixed effect for time-variant covariate
    eff.score.beta <- t(X) %*% (y - cond.expect.y)

    ## SIGMA
    eff.score.sigma2 <- t(zeta) %*% (cond.expect.y - y)/sigma2
    eff.score.sigma2 <- eff.score.sigma2 + (sum(y^2) - cond.expect.sum.y2)/(2*sigma2^2)

    c(eff.score.beta0, eff.score.beta1, eff.score.beta, eff.score.sigma2)
  }
}
