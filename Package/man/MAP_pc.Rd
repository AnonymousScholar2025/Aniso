% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCPrior.R
\name{MAP_pc}
\alias{MAP_pc}
\title{Calculates the MAP estimate for linear noisy observation of the field using PC priors on all parameters.}
\usage{
MAP_pc(
  mesh,
  lambda,
  lambda1,
  lambda_epsilon,
  lambda_u,
  y,
  A,
  m_u,
  max_iterations = 300,
  log_sigma_epsilon = NULL,
  theta0 = c(-0.5, c(0.1, 0.1), 0, -3)
)
}
\arguments{
\item{mesh}{The mesh}

\item{lambda}{A hyperparameter controlling the size of kappa.}

\item{lambda1}{A hyperparameter controlling the size of |v|.}

\item{lambda_epsilon}{A hyperparameter controlling the size of sigma_epsilon.}

\item{y}{A vector with length equal to the number of basis elements n representing the observed data.}

\item{A}{Matrix of size nxn representing the transformation A}

\item{m_u}{A vector with length n representing the prior mean m_u}

\item{max_iterations}{Maximum number of iterations for optim, by default 300}

\item{log_sigma_epsilon}{Variance of noise, if NULL, it is estimated by the MAP}
}
\value{
The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
}
\description{
Calculated by maximizing log posterior using optim. Only stationary parameters are accepted.
}
