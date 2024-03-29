% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cox_simulation.R
\name{cox_simulation}
\alias{cox_simulation}
\title{Simulate survival response based on Cox PH model}
\usage{
cox_simulation(
  n,
  logHR,
  A,
  dist = "Exponential",
  lambda = 1,
  nu = NULL,
  alpha = NULL
)
}
\arguments{
\item{n}{sample size}

\item{logHR}{coefficients for linear predictor of Cox model}

\item{A}{matrix of auxiliary covariates for linear predictor of Cox model}

\item{dist}{desired outcome distribution; choices include "Exponential," (default) "Weibull," or "Gompertz"}

\item{lambda}{scale parameter for Exponential, Weibull, and Gompertz distributions. Must be positive. Default is \code{1}}

\item{nu}{shape parameter for Weibull distributions. Must be positive. Default is \code{NULL}}

\item{alpha}{shape parameter for Gompertz distribution. Default is \code{NULL}}
}
\value{
a response vector of length \code{n} generated by the specified Cox model
}
\description{
Simulate survival response based on Cox PH model that follows either
exponential, Weibull, or Gompertz distribution. Linear predictor for
the Cox model is given by \eqn{covariate * logHR}
}
