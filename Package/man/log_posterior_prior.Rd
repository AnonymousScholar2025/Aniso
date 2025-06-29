% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCPrior.R
\name{log_posterior_prior}
\alias{log_posterior_prior}
\title{Calculates  the log-posterior density for parameters (log_kappa,v, log_sigma_u, log_sigma_epsilon) with a general prior.}
\usage{
log_posterior_prior(
  log_prior,
  mesh,
  log_kappa,
  v,
  log_sigma_u = 0,
  log_sigma_epsilon,
  y,
  A,
  m_u
)
}
\arguments{
\item{log_prior}{A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))}

\item{mesh}{The mesh}

\item{log_kappa}{Logarithm of inverse correlation range}

\item{v}{2D vector that controls anisotropy}

\item{log_sigma_u}{Variance of field u. If unspecified, it is assumed to be 0.}

\item{y}{A vector with length m equal to the number of columns of A representing the observed data.}

\item{A}{Matrix of size mxn where n is the number of basis elements representing the transformation A}

\item{m_u}{A vector with length n representing the prior mean m_u. If a number is given, it is transformed into (m_u, m_u,..., m_u)}

\item{Q_epsilon}{A sparse matrix of size nxn representing the noise precision matrix}
}
\value{
The calculated log-posterior
}
\description{
Calculates  the log-posterior density of parameters (log(kappa),v, log(sigma_u), log(epsilon)))
given a linear noisy observation y= A*u + epsilon
Uses based on the prior density and the likelihood.
Only stationary parameters are accepted.
Value is up to an additive constant depending only on y
}
