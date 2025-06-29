# DESCRIPTION: This file contains functions to calculate the hyperparameters and priors associated with the priors
# used in the simulation study (PC, EG, uniform and beta priors) and to simulate from these priors.

#' @title Hyperparameter lambda1 for given quantile of the PC prior of |v|
#' @description Calculates the hyperparameter lambda1 for a given quantile alpha_a and large anisotropy a0 of the PC prior of |v|
#'
#' @param alpha_a A quantile in [0,1] of the PC prior of |v|
#' @param a0 A surprisingly high ratio of anisotropy
#' @export
#' @examples
#' alpha_a <- 0.01
#' a0 <- 10
#' lambda1 <- lambda1_quantile(a0 = a0, alpha_a = alpha_a)
lambda1_quantile <- function(a0, alpha_a = 0.01) {
  check_range(alpha_a, "alpha_a", 0, 1)
  check_range(a0, "a0", 1, Inf)
  r0 <- log(a0)
  lambda1 <- -log(alpha_a) / (fdist(r0) - fdist(0))
  return(lambda1)
}
#' @title Hyperparameter lambda for given quantile of the PC prior of kappa
#' @description Calculates the hyperparameter lambda for a given quantile alpha, tolerance kappa0, and lambda1 of the PC prior of kappa
#'
#' @param alpha A quantile in (0,1) of the PC prior of kappa
#' @param rho0 A surprisingly large correlation range
#' @export
#' @examples
#' alpha <- 0.01
#' rho0 <- 0.1
#' lambda1 <- 1
#' lambda <- lambda_quantile(alpha = alpha, rho0 = rho0, lambda1 = lambda1)
lambda_quantile <- function(rho0, lambda1, alpha = 0.01) {
  check_range(alpha, "alpha", 0, 1)
  kappa0 <- sqrt(8) / rho0
  product <- lambda1 * fdist(0)
  lambert <- lamW::lambertW0(product * exp(product) / alpha)
  lambda <- (lambert / fdist(0) - lambda1) / kappa0
  return(lambda)
  # Error in code get very small lambda need to fix later
  # return(5)
}

#' @title Hyperparameter for variance of gaussian field given quantiles
#' @description Calculates  the hyperparameter lambda such that if sigma~Exp(lambda) then sigma>sigma0 with probability alpha
#'
#' @param alpha_sigma A quantile in (0,1). By default 0.01
#' @param sigma0 A surprisingly high number
#'
#' @return The calculated hyperparameter lambda
#' @export
#' @examples
#' alpha_sigma <- 0.01
#' sigma0 <- 10
#' lambda_sigma <- lambda_variance_quantile(alpha_sigma = alpha_sigma, sigma0 = sigma0)
lambda_variance_quantile <- function(sigma0, alpha_sigma = 0.01) {
  # This warns the user if alpha is not in (0,1)
  check_range(alpha_sigma, "alpha_sigma", 0, 1)
  check_range(sigma0, "sigma0", 0, Inf)
  lambda_sigma <- -log(alpha_sigma) / sigma0
  return(lambda_sigma)
}
#' @title Parameter sigma^2 such that if v1,v2 ~ N(0, sigma^2) the anisotropy ratio exp(|v|) is larger than a0 with probability alpha
#' @description Calculates the parameter sigma^2 such that if v1,v2 are iid N(0, sigma^2), then the anisotropy ratio exp(|v|) is larger than a0 with probability alpha
#' Uses that |v|^2 ~ Exp(1/2sigma^2)
#'
#' @param alpha_v A quantile in (0,1). By default 0.01
#' @param a0 A surprisingly high ratio of anisotropy
#'
#' @return The calculated parameter sigma^2
#' @export
#' @examples
#' alpha_v <- 0.01
#' a0 <- 10
#' sigma2_v <- sigma2_quantile_v(alpha_v = alpha_v, a0 = a0)
sigma2_quantile_v <- function(a0, alpha_v = 0.01) {
  check_range(alpha_v, "alpha_v", 0, 1)
  check_range(a0, "a0", 1, Inf)
  r02 <- log(a0)^2
  # r^2 follows a exponential distribution with rate 1/2sigma^2
  sigma2_v <- -r02 / (2 * log(alpha_v))
  return(sigma2_v)
}
#' @title Parameter lambda_k such that if kappa~Exp(lambda_k) the correlation range sqrt{8}/kappa is smaller than rho0 with probability alpha
#' @description Calculates  the parameter lambda_k such that if kappa~Exp(lambda_k) then the correlation range sqrt{8}/kappa is smaller than rho0 with probability alpha
#'
#' @param alpha_k A quantile in (0,1), by default 0.01
#' @param rho0 A surprisingly small correlation range
#'
#' @return The calculated parameter lambda_k
#' @export
#' @examples
#' alpha_k <- 0.01
#' rho0 <- 0.1
#' lambda_k <- lambda_quantile_kappa(alpha_k = alpha_k, rho0 = rho0)
lambda_quantile_kappa <- function(rho0, alpha_k = 0.01) {
  check_range(alpha_k, "alpha_k", 0, 1)
  check_range(rho0, "rho0", 0, Inf)
  kappa0 <- sqrt(8) / rho0
  lambda_k <- -log(alpha_k) / kappa0
  return(lambda_k)
}

#' @title Log prior on anisotropy, noise and variance of field supposing v1,v2~N(0,sigma_v^2) and kappa~Exp(lambda_k)
#' and with PC priors on noise and variance of field, given certain quantiles.
#'
#' @description
#' Calculates  the log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' supposing v1,v2~N(0,sigma_v^2) and kappa~Exp(lambda_k)
#' and with PC priors on noise and variance of field, given certain quantiles.
#' such that the anisotropy ratio exp(|v|) is larger than a0 with probability alpha and the correlation range sqrt{8}/kappa is smaller than rho0 with probability alpha
#'
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param alpha_k A quantile in (0,1) for kappa. If NULL, alpha_k=alpha
#' @param alpha_v A quantile in (0,1) for anisotropy |v|. If NULL, alpha_v=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#'
#' @return The calculated log of the prior on  theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' alpha_k <- 0.01
#' alpha_v <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' log_prior_gaussian <- log_EG_prior_quantile(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, alpha_k = alpha_k, alpha_v = alpha_v, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0)
log_EG_prior_quantile <- function(sigma_u0, sigma_epsilon0, a0, rho0, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, alpha_k = NULL, alpha_v = NULL) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  if (is.null(alpha_k)) {
    alpha_k <- alpha
  }
  if (is.null(alpha_v)) {
    alpha_v <- alpha
  }
  # Calculate the hyperparameters of the priors
  sigma2_v <- sigma2_quantile_v(alpha_v = alpha_v, a0 = a0)
  lambda_k <- lambda_quantile_kappa(alpha_k = alpha_k, rho0 = rho0)
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)

  # Defines the log prior
  log_prior <- function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
    # Logarithm of log exponential prior on kappa
    log_kappa_term <- log(lambda_k) - lambda_k * exp(log_kappa) + log_kappa
    v_term <- log_gaussian_density(x = v[1], mu = 0, log_sigma = log(sigma2_v) / 2) + log_gaussian_density(x = v[2], mu = 0, log_sigma = log(sigma2_v) / 2) # log_gaussian takes standard deviation not variance
    log_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
    log_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
    log_kappa_term + v_term + log_noise_value + log_sigma_u_value
  }
  return(log_prior)
}
#' @title Exponential gaussian prior on absolute value of v
#' @description Calculates  the prior on the absolute value of |v|, supposing vi~N(0,sigma^2)
#' @param alpha_v A quantile in (0,1) for the absolute value of v
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @return A function for the the prior on the absolute value of v
#' @export

EG_prior_r <- function(alpha_v, a0) {
  check_range(alpha_v, "alpha_v", 0, 1)
  check_range(a0, "a0", 1, Inf)
  sigma2_v <- sigma2_quantile_v(alpha_v = alpha_v, a0 = a0)
  # The absolute value follows a rayleigh distribution with parameter sigma2_v
  prior <- function(r) {
    r / sigma2_v * exp(-r^2 / (2 * sigma2_v))
  }
  return(prior)
}

#' @title Log uniform density
#' @description Calculates the log of the uniform density on [a,b]
#' @param a A number
#' @param b A number
#' @return The function for the log of the uniform density on [a,b]
#' @export
#' @examples
#' a <- 0
#' b <- 1
#' log_uniform_density <- log_uniform_density(a = a, b = b)
#' log_uniform_density(0.5)
log_uniform_density <- function(a, b) {
  log_uniform_density <- function(x) {
    if (x < a || x > b) {
      return(-Inf)
    } else {
      return(-log(b - a))
    }
  }
  return(log_uniform_density)
}
#' @title Calculate support of the uniform prior
#' @description Calculates the support of the uniform prior for theta.
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#' @param L A surprisingly large correlation range (e.g. size of domain), in (0,infinity)
#' @param width_support_factor A factor to multiply the width of the support of the uniform distribution. By default 2
#' @return A list of upper and lower bounds for the support of the uniform prior. Each sublist is of length 5, corresponding to log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon
#' @export
#' @examples
#' a0 <- 10
#' rho0 <- 0.1
#' L <- 100
#' width_support_factor <- 2
#' support <- support_uniform(a0 = a0, rho0 = rho0, L = L, width_support_factor = width_support_factor)
support_uniform <- function(a0, rho0, L, width_support_factor = 2) {
  kappa0 <- sqrt(8) / rho0
  log_kappa_support <- c(-log(width_support_factor * L), -log(rho0 / width_support_factor)) + 1 / 2 * log(8)
  abs_v_support <- c(0, log(width_support_factor * a0) / sqrt(2))
  log_sigma_u_support <- c(-Inf, Inf)
  log_sigma_epsilon_support <- c(-Inf, Inf)
  lower <- c(log_kappa_support[1], -abs_v_support[1], -abs_v_support[1], log_sigma_u_support[1], log_sigma_epsilon_support[1])
  upper <- c(log_kappa_support[2], abs_v_support[2], abs_v_support[2], log_sigma_u_support[2], log_sigma_epsilon_support[2])
  return(list(lower = lower, upper = upper))
}

#' @title Log uniform prior on anisotropy, noise and variance of field
#' @description Calculates  the log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon)) supposing:
#' kappa ~ Uniform(-log(w*L)+1/2*log(8),-log(rho0/w)+1/2*log(8)), so that rho= sqrt{8}/kappa has support [rho0/w,w*L]
#'  |v1|,|v2| ~ Uniform(0,log(w*a0)/sqrt(2)), so that exp(|v|) has support [1,a0*w]
#' and with PC priors on noise and variance of field, given certain quantiles.
#'
#' @param rho0 A surprisingly small correlation range, correlation is not allowed to go below rho0/2
#' @param L A surprisingly large correlation range (e.g. size of domain), correlation is not allowed to go above L
#' @param width_support_factor A factor to multiply the width of the support of the uniform distribution. By default 2. If Inf,an improper uniform prior is used
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity), anisotropy ratio is not allowed to go below a0/2 or above 2a0.
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#'
#' @return The calculated log of the prior on  theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' L <- 100
#'
#' log_prior <- log_prior_uniform(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, width_support_factor = 2)
log_prior_uniform <- function(sigma_u0, sigma_epsilon0, a0 = NULL, rho0 = NULL, L = NULL, width_support_factor = 2, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  # This warns the user if alpha is not in (0,1)
  if (alpha <= 0 || alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  if (alpha_u <= 0 || alpha_u >= 1) {
    warning("alpha_u should be in (0,1)")
  }
  if (alpha_epsilon <= 0 || alpha_epsilon >= 1) {
    warning("alpha_epsilon should be in (0,1)")
  }

  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)
  check_range(L, "L", 0, Inf)

  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)

  # Defines the log prior
  log_prior <- function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
    # Terms for variance
    lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
    lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)
    variance_term <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon) + log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
    if (width_support_factor == Inf) {
      return(variance_term)
    }
    # Terms for anisotropy
    log_uniform_density_log_kappa <- log_uniform_density(a = -log(width_support_factor * L) + 1 / 2 * log(8), b = -log(rho0 / width_support_factor) + 1 / 2 * log(8))
    log_uniform_density_v <- log_uniform_density(a = 0, b = log(width_support_factor * a0) / sqrt(2))
    log_kappa_term <- log_uniform_density_log_kappa(log_kappa)
    v_term <- 0.5 * (log_uniform_density_v(abs(v[1])) + log_uniform_density_v(abs(v[2])))


    return(log_kappa_term + v_term + variance_term)
  }
  return(log_prior)
}
#' @title Log beta density
#' @description Calculates the log of the beta density on [a,b] centred at the midpoint with shape parameters shape
#' @param a A number
#' @param b A number
#' @param shape A shape parameter of the beta distribution
#' @return The function for the log of the beta density on [a,b]
#' @export
#' @examples
#' a <- 0
#' b <- 1
#' shape <- 1.1
#' log_beta_density <- log_beta_density(a = a, b = b, shape = shape)
log_beta_density <- function(a, b, shape) {
  log_beta_density <- function(x) {
    x_scaled <- (x - a) / (b - a)
    log_density <- log(dbeta(x_scaled, shape1 = shape, shape2 = shape))
    return(log_density - log(b - a))
  }
  return(log_beta_density)
}
#' @title Log beta prior on anisotropy, noise and variance of field
#' @description Calculates  the log of the prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon)) supposing:
#' log_kappa ~ Beta(shape,shape) such that rho=sqrt{8}/kappa has support [rho0/w,w*L]
#' v1,v2 ~ Beta(shape,shape) such that exp(|v|) has support [1,a0*w]
#' Both shape parameters of the beta distribution are set to shape
#'
#' @param rho0 A surprisingly small correlation range, correlation is not allowed to go below rho0/2
#' @param L A surprisingly large correlation range (e.g. size of domain), correlation is not allowed to go above L
#' @param width_support_factor A factor to multiply the width of the support of the uniform distribution. By default 2
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity), anisotropy ratio is not allowed to go below a0/2 or above 2a0
#' @param shape A shape parameter of the beta distribution
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#'
#' @return The calculated log of the prior on  theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' L <- 100
#' width_support_factor <- 2
#' shape <- 1.1
#'
#' log_prior <- log_prior_beta(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, shape = shape, width_support_factor = width_support_factor)
log_prior_beta <- function(sigma_u0, sigma_epsilon0, a0, rho0, L, shape, width_support_factor = 2, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  # This warns the user if alpha is not in (0,1)
  if (alpha <= 0 || alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  if (alpha_u <= 0 || alpha_u >= 1) {
    warning("alpha_u should be in (0,1)")
  }
  if (alpha_epsilon <= 0 || alpha_epsilon >= 1) {
    warning("alpha_epsilon should be in (0,1)")
  }

  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)
  check_range(L, "L", 0, Inf)

  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)

  # Defines the log prior
  log_prior <- function(log_kappa, v, log_sigma_u, log_sigma_epsilon) {
    # Terms for anisotropy
    log_beta_density_log_kappa <- log_beta_density(a = -log(width_support_factor * L) + 1 / 2 * log(8), b = -log(rho0 / width_support_factor) + 1 / 2 * log(8), shape = shape)
    support_v <- log(width_support_factor * a0) / sqrt(2)
    log_beta_density_v <- log_beta_density(a = -support_v, b = support_v, shape = shape)
    log_kappa_term <- log_beta_density_log_kappa(log_kappa)
    v_term <- log_beta_density_v(v[1]) + log_beta_density_v(v[2])
    log_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
    log_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
    log_kappa_term + v_term + log_noise_value + log_sigma_u_value
  }
  return(log_prior)
}


#' @title Log PC prior on theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon)) given certain quantiles.
#' @description Calculates  the log of the PC prior on the (log(kappa),v, log(sigma_u), log(sigma_epsilon)) given certain quantiles.
#'
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise
#' @param alpha_k A quantile in (0,1) for kappa
#' @param alpha_v A quantile in (0,1) for anisotropy |v|
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#'
#' @return The calculated log of the PC prior on  theta=(log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' @export
#' @examples
#' alpha <- 0.05
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' log_pc_prior <- log_pc_prior_quantile(sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, alpha = alpha)
log_pc_prior_quantile <- function(sigma_u0, sigma_epsilon0, a0, rho0, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, alpha_k = NULL, alpha_v = NULL) {
  # This warns the user if alpha is not in (0,1)
  if (alpha <= 0 || alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  if (is.null(alpha_k)) {
    alpha_k <- alpha
  }
  if (is.null(alpha_v)) {
    alpha_v <- alpha
  }
  check_range(alpha, "alpha", 0, 1)
  check_range(alpha_u, "alpha_u", 0, 1)
  check_range(alpha_epsilon, "alpha_epsilon", 0, 1)
  check_range(alpha_k, "alpha_k", 0, 1)
  check_range(alpha_v, "alpha_v", 0, 1)
  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)

  # Calculate the hyperparameters of the PC priors for anisotropy
  lambda1 <- lambda1_quantile(a0 = a0, alpha_a = alpha_v)
  lambda <- lambda_quantile(rho0 = rho0, lambda1 = lambda1, alpha = alpha_k)
  # Calculate the hyperparameters of the PC priors for variances
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)

  # Defines the log prior
  log_pc_prior <-
    function(log_kappa,
             v,
             log_sigma_u,
             log_sigma_epsilon) {
      log_pc_aniso_value <-
        log_pc_prior_aniso(
          lambda = lambda,
          lambda1 = lambda1,
          log_kappa = log_kappa,
          v = v
        )
      log_pc_noise_value <-
        log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)
      log_pc_sigma_u_value <-
        log_pc_prior_noise_variance(lambda_epsilon = lambda_u, log_sigma_epsilon = log_sigma_u)
      log_pc_value <-
        log_pc_aniso_value + log_pc_noise_value + log_pc_sigma_u_value
      return(log_pc_value)
    }
  return(log_pc_prior)
}

#' @title Log PC prior on r using quantiles
#' @description Calculates  the log of the PC prior on r given quantiles
#'
#' @param alpha A quantile in (0,1) for anisotropy
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @return The PC prior function on r
#' @export

PC_prior_r_quantile <- function(a0, alpha) {
  check_range(alpha, "alpha_r", 0, 1)
  check_range(a0, "a0", 1, Inf)
  lambda1 <- lambda1_quantile(a0 = a0, alpha_a = alpha)
  prior <- function(r) {
    PC_prior_r(r = r, lambda1 = lambda1)
  }
}



#' @title Helper function to check range of parameters
#' @description Checks if the parameters are in the correct range
#' @param x A number
#' @param name A string
#' @param lower A number
#' @param upper A number
#' @return A warning if x is not in the range [lower,upper]
#' @export
#' @examples
#' check_range(x = 0.5, name = "x", lower = 0, upper = 1)
check_range <- function(x, name, lower, upper) {
  if (x < lower || x > upper) {
    warning(paste(name, "should be in (", lower, ",", upper, ")"))
  }
}

#' @title Simulation of theta given quantiles of the PC prior
#' @description Simulates theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon)) from the PC prior given quantiles of the PC prior
#'
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise
#' @param alpha_k A quantile in (0,1) for kappa
#' @param alpha_v A quantile in (0,1) for anisotropy |v|
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#' @param m Number of samples, by default 1.
#'
#' @return A list with four elements: log_kappa, v, log_sigma_u, log_sigma_epsilon
#' @export
#' @examples
#' alpha <- 0.05
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' m <- 10
#' result <- sim_theta_pc_quantile(alpha = alpha, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = m)
sim_theta_pc_quantile <- function(sigma_u0, sigma_epsilon0, a0, rho0, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, alpha_k = NULL, alpha_v = NULL, m = 1) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  if (is.null(alpha_k)) {
    alpha_k <- alpha
  }
  if (is.null(alpha_v)) {
    alpha_v <- alpha
  }
  check_range(alpha, "alpha", 0, 1)
  check_range(alpha_u, "alpha_u", 0, 1)
  check_range(alpha_epsilon, "alpha_epsilon", 0, 1)
  check_range(alpha_k, "alpha_k", 0, 1)
  check_range(alpha_v, "alpha_v", 0, 1)
  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)

  # Calculate the hyperparameters of the PC priors for anisotropy
  lambda1 <- lambda1_quantile(a0 = a0, alpha_a = alpha_v)
  lambda <- lambda_quantile(rho0 = rho0, lambda1 = lambda1, alpha = alpha_k)
  # Calculate the hyperparameters of the PC priors for variances
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)
  return(sim_theta_pc(lambda = lambda, lambda1 = lambda1, lambda_u = lambda_u, lambda_epsilon = lambda_epsilon, m = m))
}

#' @title Simulation of theta given quantiles for the uniform prior
#' @description Simulates theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon)) from the uniform prior given quantiles for the uniform prior
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#' @param L A surprisingly large correlation range (e.g. size of domain), in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @return A list with four elements: log_kappa, v, log_sigma_u, log_sigma_epsilon
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' L <- 100
#' m <- 10
#' result <- sim_theta_uniform(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L)
sim_theta_uniform <- function(sigma_u0, sigma_epsilon0, a0, rho0, L, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, width_support_factor = 2) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  # This warns the user if alpha is not in (0,1)
  if (alpha <= 0 || alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  if (alpha_u <= 0 || alpha_u >= 1) {
    warning("alpha_u should be in (0,1)")
  }
  if (alpha_epsilon <= 0 || alpha_epsilon >= 1) {
    warning("alpha_epsilon should be in (0,1)")
  }
  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)
  check_range(L, "L", 0, Inf)

  log_kappa <- runif(n = 1, min = -log(width_support_factor * L) + 1 / 2 * log(8), max = -log(rho0 / width_support_factor) + 1 / 2 * log(8))
  abs_v1 <- runif(n = 1, min = 0, max = log(2 * a0) / sqrt(2))
  abs_v2 <- runif(n = 1, min = 0, max = log(2 * a0) / sqrt(2))
  v1 <- sign(runif(n = 1, min = -1, max = 1)) * abs_v1
  v2 <- sign(runif(n = 1, min = -1, max = 1)) * abs_v2
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)
  sigma_u <- rexp(1, rate = lambda_u)
  sigma_epsilon <- rexp(1, rate = lambda_epsilon)
  return(list(log_kappa = log_kappa, v = cbind(v1, v2), log_sigma_u = log(sigma_u), log_sigma_epsilon = log(sigma_epsilon)))
}

#' @title Simulation of theta given quantiles for the beta prior
#' @description Simulates theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon)) from the beta prior. That is:
#' correlation range rho = sqrt(8)/kappa ~ beta(shape,shape) with support (rho0/2,L),
#' v1,v2 ~beta(shape,shape)with support (log(a0/w)/sqrt(2),log(wa0)/sqrt(2)), and with PC priors on noise and variance of field, given certain quantiles.
#' Both shape parameters of the beta distribution are set to shape
#'
#' @param rho0 A surprisingly small correlation range, correlation is not allowed to go below rho0/2
#' @param L A surprisingly large correlation range (e.g. size of domain), correlation is not allowed to go above L
#' @param width_support_factor A factor to multiply the width of the support of the uniform distribution. By default 2
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity), anisotropy ratio is not allowed to go below a0/2 or above 2a0
#' @param shape A shape parameter of the beta distribution
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param m Number of samples, by default 1.
#' @return A list with four elements: log_kappa, v, log_sigma_u, log_sigma_epsilon
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' L <- 100
#' width_support_factor <- 2
#' shape <- 1.1
#' m <- 10
#' result <- sim_theta_beta(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, L = L, shape = shape, width_support_factor = width_support_factor, m = m)
sim_theta_beta <- function(sigma_u0, sigma_epsilon0, a0, rho0, L, shape, width_support_factor = 2, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, m = 1) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  # This warns the user if alpha is not in (0,1)
  if (alpha <= 0 || alpha >= 1) {
    warning("alpha should be in (0,1)")
  }
  if (alpha_u <= 0 || alpha_u >= 1) {
    warning("alpha_u should be in (0,1)")
  }
  if (alpha_epsilon <= 0 || alpha_epsilon >= 1) {
    warning("alpha_epsilon should be in (0,1)")
  }
  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)
  check_range(L, "L", 0, Inf)

  log_kappa <- rbeta(m, shape, shape) * (log(width_support_factor * L) - log(rho0 / width_support_factor)) - log(L * width_support_factor) + 1 / 2 * log(8)
  abs_v1 <- rbeta(m, shape, shape) * (log(width_support_factor * a0) / sqrt(2))
  abs_v2 <- rbeta(m, shape, shape) * (log(width_support_factor * a0) / sqrt(2))
  v1 <- sign(runif(m, min = -1, max = 1)) * abs_v1
  v2 <- sign(runif(m, min = -1, max = 1)) * abs_v2
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)
  sigma_u <- rexp(1, rate = lambda_u)
  sigma_epsilon <- rexp(1, rate = lambda_epsilon)
  return(list(log_kappa = log_kappa, v = cbind(v1, v2), log_sigma_u = log(sigma_u), log_sigma_epsilon = log(sigma_epsilon)))
}
#' @title Sim not_pc: log(kappa)~N(1,sigma_k^2) and v1,v2~N(0,sigma_v^2)
#' and with PC priors on noise and variance of field, given certain quantiles.
#'
#' @description
#' CSimulates from (log(kappa),v, log(sigma_u), log(sigma_epsilon))
#' supposing log(|v|)~N(1,sigma_v^2) and kappa~Exp(lambda_k)
#' and with PC priors on noise and variance of field, given certain quantiles.
#' such that the anisotropy ratio exp(|v|) is larger than a0 with probability alpha and the correlation range sqrt{8}/kappa is smaller than rho0 with probability alpha
#'
#' @param alpha A quantile in (0,1) for all the parameters. By default, alpha = 0.01
#' @param alpha_u A quantile in (0,1) for the variance of the field. If NULL, alpha_u=alpha
#' @param alpha_epsilon A quantile in (0,1) for the variance of the noise. If NULL, alpha_epsilon=alpha
#' @param alpha_k A quantile in (0,1) for kappa. If NULL, alpha_k=alpha
#' @param alpha_v A quantile in (0,1) for anisotropy |v|. If NULL, alpha_v=alpha
#' @param sigma_u0 A surprisingly high variance of field, in (0,infinity)
#' @param sigma_epsilon0 A surprisingly high variance of noise, in (0,infinity)
#' @param a0 A surprisingly high ratio of anisotropy, in (1,infinity)
#' @param rho0 A surprisingly small correlation range, in (0,infinity)
#' @param m Number of samples, by default 1.
#'
#' @return The simulated values of log(kappa), v, log(sigma_u), log(sigma_epsilon)
#' @export
#' @examples
#' alpha_u <- 0.01
#' alpha_epsilon <- 0.01
#' alpha_k <- 0.01
#' alpha_v <- 0.01
#' sigma_u0 <- 10
#' sigma_epsilon0 <- 1
#' a0 <- 10
#' rho0 <- 0.1
#' m <- 10
#' result <- sim_not_pc(alpha_u = alpha_u, alpha_epsilon = alpha_epsilon, alpha_k = alpha_k, alpha_v = alpha_v, sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = m)
sim_not_pc <- function(sigma_u0, sigma_epsilon0, a0, rho0, alpha = 0.01, alpha_u = NULL, alpha_epsilon = NULL, alpha_k = NULL, alpha_v = NULL, m = 1) {
  # This sets the NULL values to alpha
  if (is.null(alpha_u)) {
    alpha_u <- alpha
  }
  if (is.null(alpha_epsilon)) {
    alpha_epsilon <- alpha
  }
  if (is.null(alpha_k)) {
    alpha_k <- alpha
  }
  if (is.null(alpha_v)) {
    alpha_v <- alpha
  }
  check_range(alpha, "alpha", 0, 1)
  check_range(alpha_u, "alpha_u", 0, 1)
  check_range(alpha_epsilon, "alpha_epsilon", 0, 1)
  check_range(alpha_k, "alpha_k", 0, 1)
  check_range(alpha_v, "alpha_v", 0, 1)
  check_range(a0, "a0", 1, Inf)
  check_range(rho0, "rho0", 0, Inf)
  check_range(sigma_u0, "sigma_u0", 0, Inf)
  check_range(sigma_epsilon0, "sigma_epsilon0", 0, Inf)

  # Calculate the hyperparameters of the priors
  sigma2_v <- sigma2_quantile_v(alpha_v = alpha_v, a0 = a0)
  lambda_k <- lambda_quantile_kappa(alpha_k = alpha_k, rho0 = rho0)
  lambda_u <- lambda_variance_quantile(alpha_sigma = alpha_u, sigma0 = sigma_u0)
  lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha_epsilon, sigma0 = sigma_epsilon0)

  kappa <- rexp(m, rate = lambda_k)
  v1 <- rnorm(m, mean = 0, sd = sqrt(sigma2_v))
  v2 <- rnorm(m, mean = 0, sd = sqrt(sigma2_v))
  sigma_u <- rexp(m, rate = lambda_u)
  sigma_epsilon <- rexp(m, rate = lambda_epsilon)
  list(log_kappa = log(kappa), v = cbind(v1, v2), log_sigma_u = log(sigma_u), log_sigma_epsilon = log(sigma_epsilon))
}
