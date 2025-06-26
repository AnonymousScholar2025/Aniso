# DESCRIPTION: This script contains the functions to calculate the scores used to validate the different models
# used in the precipitation study. This includes calculation of the priors, posteriors and the leave out one square error, continuous ranked probability score and Dawid-Sebastiani score.

#' @title Log-posterior density for precipitation model
#'
#' @description
#' Calculates  the log-posterior density of parameters (log(kappa),v, log(sigma_u), log(epsilon)))
#' given a linear noisy observation y= A*u + beta_0 +beta_1 h + epsilon
#' Only stationary parameters are accepted.
#' Value is up to an additive constant depending only on y
#' beta =(beta_0,beta_1) is the linear effect and is modelled to be Gaussian independent from the rest of the parameters
#'
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param log_kappa Logarithm of inverse correlation range
#' @param v 2D vector that controls anisotropy
#' @param log_sigma_u Variance of field u. If unspecified, it is assumed to be 0.
#' @param log_sigma_epsilon Variance of noise
#' @param y A vector with length equal to the number of columns of A representing the observed data.
#' @param A Matrix of size mxn where n is the number of basis elements representing the transformation A
#' @param h A numeric vector representing the covariate
#' @param m_u A vector with length n representing the prior mean m_u. If a number is given, it is transformed into (m_u, m_u,..., m_u). By default, set to 0
#' @param m_beta A vector with length 2 representing the prior mean m_beta. If a number is given, it is transformed into (m_beta, m_beta). By default, set to 0
#' @return The calculated log-posterior

#' @return The calculated log-posterior
#' @export
log_posterior_prior_precipitation <- function(log_prior_theta, log_sigma_beta, mesh, log_kappa, v = NULL, log_sigma_u = 0, log_sigma_epsilon, y, A, h, m_u = 0, m_beta = 0) {
    kappa <- exp(log_kappa)
    sigma_epsilon <- exp(log_sigma_epsilon)
    tau_beta <- 1 / exp(log_sigma_beta)^2

    # Calculates log-prior density of theta
    ANISO <- length(formals(log_prior_theta)) == 4
    ISO <- length(formals(log_prior_theta)) == 3
    if (ANISO) {
        log_prior_value <- log_prior_theta(log_kappa, v, log_sigma_u, log_sigma_epsilon)
    } else if (ISO) {
        log_prior_value <- log_prior_theta(log_kappa, log_sigma_u, log_sigma_epsilon)
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }


    # Calculates anisotropy
    n <- nrow(mesh$loc)
    m <- length(y)
    if (m != nrow(A)) {
        stop("y and A must have the same number of rows")
    }
    if (m != length(h)) {
        stop("y and h must have the same length")
    }
    if (ISO) {
        v <- c(0, 0)
    }

    kappa_values <- rep(kappa, n)
    vec_values <- matrix(v, n, 2, byrow = TRUE)
    aniso <- list(kappa = kappa_values, vec = vec_values)

    # Calculates log-density of the distribution of u at m_ub knowing (kappa, v)

    Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
    Q_bb <- diag(rep(tau_beta, 2), nrow = 2, ncol = 2)
    Q_bu <- Matrix(0, nrow = 2, ncol = n)
    Q_b_and_u <- rbind(cbind(Q_bb, Q_bu), cbind(t(Q_bu), Q_u))

    # Calculate A matrix
    A_beta <- cbind(cbind(rep(1, m), h), A)

    # Mean of (beta,u)
    if (length(m_beta) == 1) {
        m_beta <- rep(m_beta, 2)
    }
    if (length(m_u) == 1) {
        m_u <- rep(m_u, n)
    }
    m_bu <- c(m_beta, m_u)
    ub <- m_bu

    # Contribution of field knowing parameter theta
    logGdty_prior <- logGdensity(x = ub, mu = m_bu, Q = Q_b_and_u)

    # Calculates Q_epsilon, and conditional terms Q_{ub|y,theta} and m_{ub|y,theta}
    Q_epsilon <- Matrix::Diagonal(length(y), 1 / sigma_epsilon^2)
    Q_uby_theta <- Q_b_and_u + t(A_beta) %*% Q_epsilon %*% A_beta
    m_uby_theta <- solve(Q_uby_theta, Q_b_and_u %*% m_bu + t(A_beta) %*% Q_epsilon %*% y)

    # Calculates  log-density of the posterior distribution of ub given y and theta
    logGdty_posterior <- logGdensity(x = ub, mu = m_uby_theta, Q = Q_uby_theta)

    # Calculates  log-density of the observation of y given ub, theta
    logGdty_observation <- logGdensity(x = y, mu = A_beta %*% ub, Q = Q_epsilon)

    # Calculates  log-posterior density
    unname(log_prior_value + logGdty_prior + logGdty_observation - logGdty_posterior)
}


#' @title Calculates the MAP estimate for linear noisy observation of field with a general prior on theta = (log(kappa),v, log(sigma_u), log(sigma_epsilon)).
#'
#' @description
#' Calculated by maximizing log posterior using optim. Only stationary parameters are accepted.
#' First uses Nelder-Mead to find a good starting point for the optimization
#' and then uses BFGS to maximize the log-posterior density
#'
#' @param log_prior A function that calculates the log prior of (log(kappa), v, log(sigma_u), log(sigma_epsilon)).
#' If not specified, it is assumed to be the function that returns 0.
#' @param mesh The mesh
#' @param y A vector with length equal to the number of basis elements n representing the observed data.
#' @param A Matrix of size nxn representing the transformation A
#' @param h A numeric vector representing the covariate
#' @param m_ub A vector with length n representing the prior mean m_ub
#' @param max_iterations Maximum number of iterations for optim, by default 300
#' @param log_sigma_epsilon Variance of noise, if NULL, it is estimated by the MAP
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (log(0.5), 1, 2, 1, 1)
#' @param hessian If TRUE, the hessian is returned. By default, set to TRUE
#'
#'
#' @return The parameters (log_kappa, v, log_sigma_u, log_sigma_epsilon) that maximize the posterior
#' @export
MAP_spatial_effect <- function(
    log_prior_theta, log_sigma_beta, mesh, y, A, h, m_u = 0, m_beta = 0, max_iterations = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3, 0, 0), do_u_want_hessian = TRUE) {
    ANISO <- length(formals(log_prior_theta)) == 4
    ISO <- length(formals(log_prior_theta)) == 3
    if (!ANISO && !ISO) {
        stop("The prior function must have 3 or 4 arguments")
    }
    if (ANISO) {
        log_post <- function(theta) {
            log_kappa <- theta[1]
            v <- theta[2:3]
            log_sigma_u <- theta[4]
            log_sigma_epsilon <- theta[5]
            tryCatch(
                {
                    log_posterior_prior_precipitation(
                        log_prior_theta = log_prior_theta,
                        log_sigma_beta = log_sigma_beta,
                        mesh = mesh, log_kappa = log_kappa, v = v,
                        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
                        y = y, A = A, h = h, m_u = m_u, m_beta = m_beta
                    )
                },
                error = function(e) {
                    -Inf
                }
            )
        }
    } else {
        log_post <- function(theta) {
            log_kappa <- theta[1]
            log_sigma_u <- theta[2]
            log_sigma_epsilon <- theta[3]
            tryCatch(
                {
                    log_posterior_prior_precipitation(
                        log_prior_theta = log_prior_theta,
                        log_sigma_beta = log_sigma_beta,
                        mesh = mesh, log_kappa = log_kappa,
                        log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
                        y = y, A = A, h = h, m_u = m_u, m_beta = m_beta
                    )
                },
                error = function(e) {
                    -Inf
                }
            )
        }
    }
    if (ISO) {
        if (length(theta0) == 5) {
            theta0 <- c(theta0[1], theta0[4], theta0[5])
        }
    }

    map <- optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = do_u_want_hessian)
    theta0 <- map$par
    tryCatch(
        {
            return(optim(par = theta0, fn = log_post, control = list(fnscale = -1, maxit = max_iterations / 2), hessian = do_u_want_hessian, method = "BFGS"))
        },
        error = function(e) {
            warning("Error in BFGS: ", e)
            return(map)
        }
    )
}

#' @title Log-prior density for isotropic parameters
#' @description
#' Calculates the log-prior density of parameters (log(kappa), log(sigma_u), log(sigma_epsilon)) given quantiles

#' @param sigma_u0 Quantile for standard deviation of the prior of the field u
#' @param sigma_epsilon0 Quantile for standard deviation of the noise
#' @param rho0 Quantile for the correlation range
#' @param alpha Quantile by default 0.01

log_pc_prior_iso <- function(
    sigma_u0, sigma_epsilon0, rho0, alpha = 0.01) {
    kappa0 <- sqrt(8) / rho0
    log_pc_prior <- function(log_kappa, log_sigma_u, log_sigma_epsilon) {
        # PC on kappa
        change_variable <- log_kappa
        lambda_kappa <- -rho0 / sqrt(8) * log(alpha)
        kappa <- exp(log_kappa)
        log_pc_kappa_value <- log(lambda_kappa) - lambda_kappa * kappa - change_variable # kappa ~ Exp(lambda_kappa)

        # PC on noise and variance
        lambda_sigma_u <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigma_u0)
        lambda_epsilon <- lambda_variance_quantile(alpha_sigma = alpha, sigma0 = sigma_epsilon0)
        log_pc_sigma_u_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_sigma_u, log_sigma_epsilon = log_sigma_u)
        log_pc_noise_value <- log_pc_prior_noise_variance(lambda_epsilon = lambda_epsilon, log_sigma_epsilon = log_sigma_epsilon)

        log_pc_value <- log_pc_kappa_value + log_pc_sigma_u_value + log_pc_noise_value
        return(log_pc_value)
    }
    return(log_pc_prior)
}


#' @title Prior density function for rho if field is isotropic
#' @description
#' Returns the prior density function of parameters rho given quantiles

#' @param rho0 Quantile for the correlation range
#' @param alpha Quantile by default 0.01

pc_prior_rho_iso <- function(rho0, alpha = 0.01) {
    pc_rho_iso <- function(rho) {
        # Change of variables
        kappa <- sqrt(8) / rho
        change_variable <- sqrt(8) / rho^2
        # PC on rho
        lambda_kappa <- -rho0 / sqrt(8) * log(alpha)
        lambda_kappa * exp(-lambda_kappa * kappa) * change_variable # kappa ~ exp(lambda_kappa)
    }
    return(pc_rho_iso)
}
#' @title Prior density function for kappa if field is isotropic
#' @description
#' Returns the prior density function of parameters rho given quantiles

#' @param rho0 Quantile for the correlation range
#' @param alpha Quantile by default 0.01

pc_prior_kappa_iso <- function(rho0, alpha = 0.01) {
    pc_rho_iso <- function(kappa) {
        # PC on kappa
        lambda_kappa <- -rho0 / sqrt(8) * log(alpha)
        lambda_kappa * exp(-lambda_kappa * kappa) # kappa ~ exp(lambda_kappa)
    }
    return(pc_rho_iso)
}





#' @title Log-posterior density of field u given a linear noisy observation y= A*u + beta_0 +beta_1 h + epsilon
#' @description
#' Calculates the log-posterior density of field u given a linear noisy observation y= A*u + beta_0 +beta_1 h + epsilon
#' Only stationary parameters are accepted.
#' Value is up to an additive constant depending only on y
#' beta =(beta_0,beta_1) is the linear effect and is modelled to be Gaussian independent from the rest of the parameters
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param y A vector with length equal to the number of columns of A representing the observed data.
#' @param A Matrix of size mxn where n is the number of basis elements representing the transformation A
#' @param h A numeric vector representing the covariate
#' @param m_theta Arbitrary point in calculation of the log-posterior. By default, set to 0
#' @return The calculated log-posterior

log_post_u_precipitation <- function(log_prior_theta, log_sigma_beta, mesh, y, A, h, m_theta = rep(0, 5)) {
    log_post_u <- function(u) {
        kappa <- exp(log_kappa)
        sigma_epsilon <- exp(log_sigma_epsilon)
        tau_beta <- 1 / exp(log_sigma_beta)^2

        # Calculates log-prior density of theta
        log_prior_value <- log_prior_theta(log_kappa, v, log_sigma_u, log_sigma_epsilon)

        # Calculates anisotropy
        n <- nrow(mesh$loc)
        m <- length(y)
        if (m != nrow(A)) {
            stop("y and A must have the same number of rows")
        }
        if (m != length(h)) {
            stop("y and h must have the same length")
        }
        kappa_values <- rep(kappa, n)
        vec_values <- matrix(v, n, 2, byrow = TRUE)
        aniso <- list(kappa = kappa_values, vec = vec_values)

        # Calculates log-density of the distribution of u at m_ub knowing (kappa, v)

        Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
        Q_bb <- diag(rep(tau_beta, 2), nrow = 2, ncol = 2)
        Q_bu <- Matrix(0, nrow = 2, ncol = n)
        Q_b_and_u <- rbind(cbind(Q_bb, Q_bu), cbind(t(Q_bu), Q_u))

        # Calculate A matrix
        A_beta <- cbind(cbind(rep(1, m), h), A)

        # Mean of (beta,u)
        if (length(m_beta) == 1) {
            m_beta <- rep(m_beta, 2)
        }
        if (length(m_u) == 1) {
            m_u <- rep(m_u, n)
        }
        m_bu <- c(m_beta, m_u)
        ub <- m_bu

        # Contribution of field knowing parameter theta
        logGdty_prior <- logGdensity(x = ub, mu = m_bu, Q = Q_b_and_u)

        # Calculates Q_epsilon, and conditional terms Q_{ub|y,theta} and m_{ub|y,theta}
        Q_epsilon <- Matrix::Diagonal(length(y), 1 / sigma_epsilon^2)
        Q_uby_theta <- Q_b_and_u + t(A_beta) %*% Q_epsilon %*% A_beta
        m_uby_theta <- solve(Q_uby_theta, Q_b_and_u %*% m_bu + t(A_beta) %*% Q_epsilon %*% y)

        # Calculates  log-density of the posterior distribution of ub given y and theta
        logGdty_posterior <- logGdensity(x = ub, mu = m_uby_theta, Q = Q_uby_theta)

        # Calculates  log-density of the observation of y given ub, theta
        logGdty_observation <- logGdensity(x = y, mu = A_beta %*% ub, Q = Q_epsilon)

        # Calculates  log-posterior density
        unname(log_prior_value + logGdty_prior + logGdty_observation - logGdty_posterior)
    }
    log_post_u
}

#' @title Smoothed log importance weights
#' @description Calculates the smoothed self normalized weights where the importance density is a Gaussian distribution
#' @param log_posterior_density A function that calculates the log posterior density of the parameters
#' @param mu_Gaussian_median The median of the Gaussian importance density of the parameters
#' @param Q_Gaussian_median The precision matrix of the Gaussian importance density of the parameters
#' @param n_weights The number of importance weights
#' @param alpha The quantile for the credible interval, by default 0.05
#' @return The smoothed log importance weights and samples

log_smoothed_importance_weights_and_samples <- function(log_posterior_density, mu_Gaussian_median, Q_Gaussian_median, n_weights, alpha = 0.05) {
    # Importance sample
    theta_sim_importance <- MASS::mvrnorm(n = n_weights, mu = mu_Gaussian_median, Sigma = solve(Q_Gaussian_median))

    # Importance weights
    log_Gaussian_median_density <- function(theta) {
        logGdensity(
            x = theta, mu = mu_Gaussian_median, Q = Q_Gaussian_median
        )
    }

    log_ratio_function <- function(theta) {
        log_posterior_density(theta) - log_Gaussian_median_density(theta)
    }
    log_importance_ratios <- apply(theta_sim_importance, 1, log_ratio_function)

    # We remove samples and ratios with - infinite importance ratios and subtract the max to avoid numerical issues as it shouldn't change the result
    theta_sim_importance <- theta_sim_importance[!is.infinite(log_importance_ratios), ]
    log_importance_ratios <- log_importance_ratios[!is.infinite(log_importance_ratios)]
    log_importance_ratios <- log_importance_ratios - max(log_importance_ratios)

    # Smoothing
    psis_result <- psis(log_importance_ratios, r_eff = NA)
    log_weights_smoothed <- psis_result$log_weights
    return(list(log_weights_smoothed = log_weights_smoothed, theta_sim_importance = theta_sim_importance))
}

importance_confidence <- function(f_values, theta_sim_importance, normalized_weights, alpha = alpha) {
    N <- length(f_values)
    integral <- mean(f_values * normalized_weights)
    variance <- mean(weights^2 * (f_values - integral)^2) / N
    symmetric_quantile <- qnorm(1 - alpha / 2)
    confidence_interval <- integral + symmetric_quantile * sqrt(variance)
}


#' @title LOO for the precipitation model
#' @description
#' Calculates the leave out one cross validation for the precipitation model
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param Y A vector with length equal to the number of columns of A representing the observed data.
#' @param A Matrix of size mxn where n is the number of basis elements representing the transformation A
#' @param h A numeric vector representing the covariate at the observed data
#' @param h_mesh A numeric vector representing the covariate at the mesh points
#' @param max_iterations_MAP Maximum number of iterations for MAP, by default 300
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (-0.5, c(0.1, 0.1), 0, -3, 0, 0)
#' @param n_weights The number of importance weights, by default 1000
#' @param path_data Path to the data
#' @param path_prediction_txt Path to the prediction
#' param alpha The quantile for the credible interval, by default 0.05
#' @return The calculated leave out one cross validation
#' @export

LOO_CV_precipitation <- function(log_prior_theta,
                                 log_sigma_beta,
                                 mesh,
                                 y,
                                 A,
                                 h,
                                 max_iterations_MAP = 300,
                                 theta0 = c(-0.5, c(0.1, 0.1), 0, -3, 0, 0),
                                 n_weights = 1000,
                                 h_mesh,
                                 path_data = NULL,
                                 path_prediction_txt = NULL,
                                 quiet = FALSE,
                                 alpha = 0.05) {
    # Dimensions for later
    n_y <- length(y)
    n_u <- nrow(mesh$loc)
    if (n_y != nrow(A)) {
        stop("y and A must have the same number of rows")
    }
    if (n_y != length(h)) {
        stop("y and h must have the same length")
    }
    ANISO <- length(formals(log_prior_theta)) == 4
    ISO <- length(formals(log_prior_theta)) == 3
    if (ANISO) {
        d_theta <- 5
        theta_true <- theta0
    } else if (ISO) {
        d_theta <- 3
        theta_true <- c(theta0[1], theta0[4], theta0[5])
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }

    # Get Gaussian median approximation of p(theta|y)
    map_precip <- MAP_spatial_effect(
        log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
        y = y, A = A, h, m_u = 0, m_beta = 0, max_iterations = max_iterations_MAP,
        theta0 = theta0
    )
    mu_Gaussian_median <- map_precip$par
    Q_Gaussian_median <- -map_precip$hessian
    # mu_Gaussian_median <- c(10, -0.43268029, 0.02833754, -0.45326264, -1.95864792) # For testing
    # map_precip <- list(par = mu_Gaussian_median)
    # Q_Gaussian_median <- 1 * diag(5)
    credible_interval <- calculate_credible_intervals_gaussian(mu_Gaussian_median, sqrt(diag(solve(Q_Gaussian_median))), alpha)
    parameter_in_CI_vector <- parameter_within_credible_intervals(theta_true, credible_interval)
    interval_scores_parameters <- interval_score(theta_true, credible_interval, alpha = alpha)
    # credible_interval <- convert_to_range_and_non_log(t(credible_interval))
    # Importance weights and samples
    if (ANISO) {
        log_posterior <- function(theta) {
            log_posterior_prior_precipitation(
                log_prior_theta = log_prior_theta,
                log_sigma_beta = log_sigma_beta,
                mesh = mesh, log_kappa = theta[1], v = theta[2:3],
                log_sigma_u = theta[4], log_sigma_epsilon = theta[5],
                y = y, A = A, h = h, m_u = 0, m_beta = 0
            )
        }
    } else if (ISO) {
        log_posterior <- function(theta) {
            log_posterior_prior_precipitation(
                log_prior_theta = log_prior_theta,
                log_sigma_beta = log_sigma_beta,
                mesh = mesh, log_kappa = theta[1],
                log_sigma_u = theta[2], log_sigma_epsilon = theta[3],
                y = y, A = A, h = h, m_u = 0, m_beta = 0
            )
        }
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }

    log_weights_and_samples <- log_smoothed_importance_weights_and_samples(
        log_posterior_density = log_posterior, mu_Gaussian_median = mu_Gaussian_median,
        Q_Gaussian_median = Q_Gaussian_median, n_weights = n_weights
    )
    log_weights <- log_weights_and_samples$log_weights_smoothed
    weights_normalized <- normalize_log_weights(log_weights)
    n_eff <- 1 / sum(weights_normalized^2)
    theta_sim_importance <- log_weights_and_samples$theta_sim_importance

    # Lists to store values
    means_y_i__y_no_i_theta_s <- matrix(0, n_y, n_weights)
    variances_y_i__y_no_i_theta_s <- matrix(0, n_y, n_weights)
    m_ub__y_theta_s <- matrix(0, 2 + n_u, n_weights) # Mean of (beta,u) given y and theta_s
    pred_at_mesh__y2 <- matrix(0, n_u, n_weights) # predicted precipitation at mesh points given y and theta_s
    tau_beta <- 1 / exp(log_sigma_beta)^2
    pb <- progress_bar$new(
        format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
        total = length((theta_sim_importance[, 1])),
        complete = "=",
        incomplete = "-",
        current = ">",
        clear = TRUE,
        width = 100
    ) # Width of the progress bar

    for (j in seq_along(theta_sim_importance[, 1])) {
        pb$tick()
        theta_s <- theta_sim_importance[j, ]
        # Calculates variance Q_u_beta|y,theta_s
        kappa <- exp(theta_s[1])
        if (ANISO) {
            v <- theta_s[2:3]
            log_sigma_u <- theta_s[4]
            sigma_epsilon <- exp(theta_s[5])
        } else {
            v <- c(0, 0)
            log_sigma_u <- theta_s[2]
            sigma_epsilon <- exp(theta_s[3])
        }
        kappa_values <- rep(kappa, n_u)
        vec_values <- matrix(v, n_u, 2, byrow = TRUE)
        aniso <- list(kappa = kappa_values, vec = vec_values)
        Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
        Q_bb <- diag(rep(tau_beta, 2), nrow = 2, ncol = 2)
        Q_bu <- Matrix(0, nrow = 2, ncol = n_u)
        Q_b_and_u <- rbind(cbind(Q_bb, Q_bu), cbind(t(Q_bu), Q_u))
        A_beta <- cbind(cbind(rep(1, n_y), h), A)
        Q_ub__y_theta <- Q_b_and_u + t(A_beta) %*% A_beta / sigma_epsilon^2
        V <- rowSums(A_beta * (A_beta %*% inla.qinv(Q_ub__y_theta))) # Variance of y_i given y and theta_s
        m_ub__y_theta <- solve(Q_ub__y_theta, t(A_beta) %*% y) / sigma_epsilon^2
        eta <- A_beta %*% m_ub__y_theta

        # Store posterior mean of u and prediction using u for parameter theta_j
        m_ub__y_theta_s[, j] <- as.vector(m_ub__y_theta)
        pred_at_mesh__y2[, j] <- m_ub__y_theta[1] + m_ub__y_theta[2] * h_mesh + m_ub__y_theta[3:(2 + n_u)]

        # Calculate and store predictive mean and variance for y_i given y_-i and theta_s
        Sigmas_y_i__y_no_i_theta <- sigma_epsilon^2 / (1 - V / sigma_epsilon^2) # predictive variance for y_i
        means_y_i__y_no_i <- y + (eta - y) / (1 - V / sigma_epsilon^2) # predictive mean for y_i
        variances_y_i__y_no_i_theta_s[, j] <- Sigmas_y_i__y_no_i_theta
        means_y_i__y_no_i_theta_s[, j] <- as.vector(means_y_i__y_no_i)
    }
    # We calculate the list of length n_y of importance approximations to the expected value of the posterior density of y_i given y_-i using the normalized weights
    m_ub__y <- apply(m_ub__y_theta_s, 1, function(x) sum(x * weights_normalized))
    pred_at_mesh__y <- apply(pred_at_mesh__y2, 1, function(x) sum(x * weights_normalized))
    expected_values <- apply(means_y_i__y_no_i_theta_s, 1, function(x) sum(x * weights_normalized))
    sigma_squared_hat <- apply(variances_y_i__y_no_i_theta_s, 1, function(x) sum(x * weights_normalized))

    # Leave out one error
    loo_error_list_squared <- (expected_values - y)^2
    loo_error <- sqrt(mean(loo_error_list_squared))


    # CRPS
    CRPS_list <- scoringRules::crps_mixnorm(y = y, m = means_y_i__y_no_i_theta_s, s = variances_y_i__y_no_i_theta_s, w = matrix(rep(weights_normalized, times = n_y), nrow = n_y, byrow = TRUE))
    CRPS <- mean(CRPS_list)


    DSS_list <- sapply(1:n_y, function(i) {
        DSS(mu = expected_values[i], Q = 1 / sigma_squared_hat[i], y = y[i])
    })
    DSS <- mean(DSS_list)
    result <- list(
        m_ub__y = m_ub__y,
        pred_at_mesh__y = pred_at_mesh__y,
        means_y_i__y_no_i_theta_s = means_y_i__y_no_i_theta_s,
        variances_y_i__y_no_i_theta_s = variances_y_i__y_no_i_theta_s,
        y_hat = expected_values,
        sigma_squared_hat = sigma_squared_hat,
        loo_error_list_squared = loo_error_list_squared,
        CRPS_list = CRPS_list,
        DSS_list = DSS_list,
        loo_error = loo_error,
        CRPS = CRPS,
        DSS = DSS,
        n_eff = n_eff,
        weights_normalized = weights_normalized,
        theta_sim_importance = theta_sim_importance,
        MAP = map_precip,
        credible_interval = credible_interval,
        parameter_in_CI_vector = parameter_in_CI_vector,
        interval_scores_parameters = interval_scores_parameters
    )
    if (!is.null(path_data)) {
        saveRDS(result, file = path_data)
    }
    if (!is.null(path_prediction_txt)) {
        write.table(pred_at_mesh__y, file = path_prediction_txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    return(result)
}



#' @title LOO for the precipitation model loop
#' @description
#' Calculates the leave out one cross validation for the precipitation model
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param h_mesh A numeric vector representing the covariate at the mesh points
#' @param data A list with x and y representing the observed data
#' @param location_indexes_matrix A matrix. Each column are the indices of the data that is observed
#' @param max_iterations_MAP Maximum number of iterations for MAP, by default 300
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (-0.5, c(0.1, 0.1), 0, -3, 0, 0)
#' @param n_weights The number of importance weights, by default 1000
#' @param path_data Path to the data
#' @param path_prediction_txt Path to the prediction
#' @return The calculated leave out one cross validation
#' @export

LOO_CV_precipitation_loop <- function(log_prior_theta, log_sigma_beta, mesh, data, location_indexes_matrix, max_iterations_MAP = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3, 0, 0), n_weights = 1000, h_mesh, path_data = NULL, path_prediction_txt = NULL) {
    # Dimensions for later
    n_loops <- dim(location_indexes_matrix)[2]
    observation_locations_list <- lapply(1:n_loops, function(i) {
        cbind(data$x[location_indexes_matrix[, i]], data$y[location_indexes_matrix[, i]])
    })

    A_list <- lapply(1:n_loops, function(i) {
        fmesher::fm_basis(mesh, loc = observation_locations_list[[i]])
    })

    # Normalize precipitation and height to be approximately 1
    y_vector <- sapply(1:n_loops, function(i) {
        data$precipitation[location_indexes_matrix[, i]] / 1000
    })
    h_vector <- sapply(1:n_loops, function(i) {
        data$altitude[location_indexes_matrix[, i]] / 1000
    })
    n_loop <- dim(location_indexes_matrix)[2]
    n_y <- dim(location_indexes_matrix)[1]
    n_u <- nrow(mesh$loc)
    if (length(formals(log_prior_theta)) == 4) {
        d_theta <- 5
    } else if (length(formals(log_prior_theta)) == 3) {
        d_theta <- 3
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }

    # Apply LOO_CV_precipitation to each observation

    pb <- progress_bar$new(
        format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
        total = n_loops,
        complete = "=",
        incomplete = "-",
        current = ">",
        clear = TRUE,
        width = 100
    )
    results <- lapply(1:n_loops, function(i) {
        pb$tick()
        tryCatch(
            {
                LOO_CV_precipitation(
                    log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
                    y = y_vector[, i], A = A_list[[i]], h = h_vector[, i], max_iterations_MAP = max_iterations_MAP,
                    theta0 = theta0, n_weights = n_weights, h_mesh = h_mesh
                )
            },
            error = function(e) {
            }
        )
    })
    not_null_indices <- sapply(results, function(x) !is.null(x))
    results <- results[not_null_indices]
    # Get scores for each loop
    loo_error_list <- sapply(results, function(x) x$loo_error)
    CRPS_list <- sapply(results, function(x) x$CRPS)
    DSS_list <- sapply(results, function(x) x$DSS)
    # Get mean of scores
    loo_error <- mean(loo_error_list)
    CRPS <- mean(CRPS_list)
    DSS <- mean(DSS_list)
    # Get bootstrap of scores
    loo_error_bootstrap_CI <- bootstrap_q(loo_error_list, B = 1000)
    CRPS_bootstrap_CI <- bootstrap_q(CRPS_list, B = 1000)
    DSS_bootstrap_CI <- bootstrap_q(DSS_list, B = 1000)

    # Get mean of MAP
    MAP <- list(
        par = rowMeans(sapply(results, function(x) x$MAP$par)),
        hessian = matrix(rowMeans(sapply(results, function(x) x$MAP$hessian)), nrow = d_theta)
    )
    MAP_non_log_par <- rowMeans(sapply(results, function(x) convert_to_range_and_non_log(t(x$MAP$par))))
    # Get mean of credible interval
    credible_interval_list <- sapply(results, function(x) x$credible_interval)
    credible_interval <- matrix(rowMeans(credible_interval_list), ncol = 2)
    # Get mean of parameter in CI
    parameter_in_CI_vector_list <- sapply(results, function(x) x$parameter_in_CI_vector)
    parameter_in_CI_vector <- rowMeans(parameter_in_CI_vector_list)
    # Get mean of interval scores
    interval_scores_parameters_list <- sapply(results, function(x) x$interval_scores_parameters)
    interval_scores_parameters <- rowMeans(interval_scores_parameters_list)

    # Get mean prediction at mesh points and of field
    pred_at_mesh__y <- rowMeans(sapply(results, function(x) x$pred_at_mesh__y))
    m_ub__y <- rowMeans(sapply(results, function(x) x$m_ub__y))
    # Save data
    result <- list(
        loo_error = loo_error,
        CRPS = CRPS,
        DSS = DSS,
        # loo_error_list_squared = loo_error_list_squared,
        loo_error_list = loo_error_list,
        CRPS_list = CRPS_list,
        DSS_list = DSS_list,
        loo_error_bootstrap_CI = loo_error_bootstrap_CI,
        CRPS_bootstrap_CI = CRPS_bootstrap_CI,
        DSS_bootstrap_CI = DSS_bootstrap_CI,
        MAP_non_log_par = MAP_non_log_par,
        credible_interval = credible_interval,
        credible_interval_list = credible_interval_list,
        parameter_in_CI_vector = parameter_in_CI_vector,
        interval_scores_parameters = interval_scores_parameters,
        pred_at_mesh__y = pred_at_mesh__y,
        m_ub__y = m_ub__y,
        MAP = MAP
    )
    if (!is.null(path_data)) {
        saveRDS(result, file = path_data)
    }
    if (!is.null(path_prediction_txt)) {
        write.table(pred_at_mesh__y, file = path_prediction_txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    return(result)
}
#' @title Leave out multiple testing
#' @description
#' Calculates the validation for the precipitation model leaving out multiple points
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param y1 A vector with length equal to the number of columns of A1 representing the observed data.
#' @param A1 Matrix of size m1xn1 where n1 is the number of basis elements representing the transformation A1
#' @param h1 A numeric vector representing the covariate at the observed data
#' @param y2 A vector with length equal to the number of columns of A2 representing the training data.
#' @param A2 Matrix of size m2xn2 where n2 is the number of basis elements representing the transformation A2
#' @param h2 A numeric vector representing the covariate at the training data
#' @param h_mesh A numeric vector representing the covariate at the mesh points
#' @param max_iterations_MAP Maximum number of iterations for MAP, by default 300
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (-0.5, c(0.1, 0.1), 0, -3, 0, 0)
#' @param n_weights The number of importance weights, by default 1000
#' @param alpha The quantile for the credible interval, by default 0.05
#' @param path_data Path to the data
#' @param path_prediction_txt Path to the prediction
#' @return The calculated leave out one cross validation
#' @export

LOM_scores <- function(log_prior_theta, log_sigma_beta, mesh, y1, A1, h1, y2, A2, h2, max_iterations_MAP = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3, 0, 0), n_weights = 1000, h_mesh, path_data = NULL, path_prediction_txt = NULL, alpha = 0.05) {
    # Dimensions for later
    n_y1 <- length(y1)
    n_y2 <- length(y2)
    n_u <- nrow(mesh$loc)
    if (n_y1 != nrow(A1) || n_y2 != nrow(A2)) {
        stop("yi and Ai must have the same number of rows")
    }
    if (n_y1 != length(h1) || n_y2 != length(h2)) {
        stop("yi and hi must have the same length")
    }

    ANISO <- length(formals(log_prior_theta)) == 4
    ISO <- length(formals(log_prior_theta)) == 3
    if (ANISO) {
        d_theta <- 5
    } else if (ISO) {
        d_theta <- 3
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }

    # Get Gaussian median approximation of p(theta|y1)
    map_precip <- MAP_spatial_effect(
        log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
        y = y1, A = A1, h1, m_u = 0, m_beta = 0, max_iterations = max_iterations_MAP,
        theta0 = theta0
    )
    mu_Gaussian_median <- map_precip$par
    Q_Gaussian_median <- -map_precip$hessian
    credible_interval <- calculate_credible_intervals_gaussian(mu_Gaussian_median, sqrt(diag(solve(Q_Gaussian_median))), alpha)
    credible_interval <- convert_to_range_and_non_log(t(credible_interval))
    # Importance weights and samples
    if (ANISO) {
        log_posterior <- function(theta) {
            log_posterior_prior_precipitation(
                log_prior_theta = log_prior_theta,
                log_sigma_beta = log_sigma_beta,
                mesh = mesh, log_kappa = theta[1], v = theta[2:3],
                log_sigma_u = theta[4], log_sigma_epsilon = theta[5],
                y = y1, A = A1, h = h1, m_u = 0, m_beta = 0
            )
        }
    } else if (ISO) {
        log_posterior <- function(theta) {
            log_posterior_prior_precipitation(
                log_prior_theta = log_prior_theta,
                log_sigma_beta = log_sigma_beta,
                mesh = mesh, log_kappa = theta[1],
                log_sigma_u = theta[2], log_sigma_epsilon = theta[3],
                y = y1, A = A1, h = h1, m_u = 0, m_beta = 0
            )
        }
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }

    log_weights_and_samples <- log_smoothed_importance_weights_and_samples(
        log_posterior_density = log_posterior, mu_Gaussian_median = mu_Gaussian_median,
        Q_Gaussian_median = Q_Gaussian_median, n_weights = n_weights
    )
    log_weights <- log_weights_and_samples$log_weights_smoothed
    weights_normalized <- normalize_log_weights(log_weights)
    n_eff <- 1 / sum(weights_normalized^2)
    theta_sim_importance <- log_weights_and_samples$theta_sim_importance

    # Lists to store values
    means_y2__y1_theta_s <- matrix(0, n_y2, n_weights)
    Sigma_y2__y1_theta_s <- array(0, dim = c(n_y2, n_y2, n_weights))
    m_ub__y1_theta_s <- matrix(0, 2 + n_u, n_weights) # Mean of (beta,u) given y and theta_s
    pred_at_mesh__y1 <- matrix(0, n_u, n_weights) # predicted precipitation at mesh points given y and theta_s
    tau_beta <- 1 / exp(log_sigma_beta)^2
    time_start <- Sys.time()
    for (j in seq_along(theta_sim_importance[, 1])) {
        print(paste0("Current loop= ", j, " out of ", n_weights))
        time_left <- (Sys.time() - time_start) / j * (n_weights - j)
        print(paste0("Time left= ", time_left))
        theta_s <- theta_sim_importance[j, ]
        # Calculates variance Q_u_beta|y,theta_s
        kappa <- exp(theta_s[1])
        if (ANISO) {
            v <- theta_s[2:3]
            log_sigma_u <- theta_s[4]
            sigma_epsilon <- exp(theta_s[5])
        } else {
            v <- c(0, 0)
            log_sigma_u <- theta_s[2]
            sigma_epsilon <- exp(theta_s[3])
        }
        kappa_values <- rep(kappa, n_u)
        vec_values <- matrix(v, n_u, 2, byrow = TRUE)
        aniso <- list(kappa = kappa_values, vec = vec_values)
        Q_u <- fm_aniso_precision(mesh, aniso, log_sigma = log_sigma_u)
        Q_bb <- diag(rep(tau_beta, 2), nrow = 2, ncol = 2)
        Q_bu <- Matrix(0, nrow = 2, ncol = n_u)
        Q_b_and_u <- rbind(cbind(Q_bb, Q_bu), cbind(t(Q_bu), Q_u))
        A1_beta <- cbind(cbind(rep(1, n_y1), h1), A1)
        A2_beta <- cbind(cbind(rep(1, n_y2), h2), A2)
        Q_ub__y1_theta <- Q_b_and_u + t(A1_beta) %*% A1_beta / sigma_epsilon^2
        m_ub__y1_theta <- solve(Q_ub__y1_theta, t(A1_beta) %*% y1) / sigma_epsilon^2
        mean_y2__y1_theta <- A2_beta %*% m_ub__y1_theta
        Sigma_y2__y1_theta <- A2_beta %*% solve(Q_ub__y1_theta) %*% t(A2_beta) + diag(x = sigma_epsilon^2, nrow = n_y2, ncol = n_y2)

        # Store posterior mean of u and prediction using u for parameter theta_j
        m_ub__y1_theta_s[, j] <- as.vector(m_ub__y1_theta)
        pred_at_mesh__y1[, j] <- m_ub__y1_theta[1] + m_ub__y1_theta[2] * h_mesh + m_ub__y1_theta[3:(2 + n_u)]

        # Calculate and store predictive mean and variance for y_1 given y_2 and theta_s
        means_y2__y1_theta_s[, j] <- as.vector(mean_y2__y1_theta)
        Sigma_y2__y1_theta_s[, , j] <- matrix(Sigma_y2__y1_theta, nrow = n_y2, ncol = n_y2)
    }
    m_ub__y1 <- apply(m_ub__y1_theta_s, 1, function(x) sum(x * weights_normalized))
    pred_at_mesh__y1 <- apply(pred_at_mesh__y1, 1, function(x) sum(x * weights_normalized))
    m_y2__y1 <- apply(means_y2__y1_theta_s, 1, function(x) sum(x * weights_normalized))
    Sigma_weighted <- sweep(Sigma_y2__y1_theta_s, 3, weights_normalized, "*")
    Sigma_y2__y1 <- apply(Sigma_weighted, c(1, 2), sum)

    # Leave out one error
    loo_error_list_squared <- (m_y2__y1 - y2)^2
    loo_error <- sqrt(mean(loo_error_list_squared))


    # CRPS is calculated as average of CRPS of marginals of y2|y1
    marginal_variances_y_i__y1_theta_s <- apply(Sigma_y2__y1_theta_s, MARGIN = 3, FUN = diag) # predictive variance for y_i
    CRPS_list <-
        scoringRules::crps_mixnorm(
            y = y2,
            m = means_y2__y1_theta_s,
            s = marginal_variances_y_i__y1_theta_s,
            w = matrix(
                rep(weights_normalized, times = n_y2),
                nrow = n_y2,
                byrow = TRUE
            )
        )
    CRPS <- mean(CRPS_list)
    DSS <- DSS(mu = m_y2__y1, Q = solve(Sigma_y2__y1), y = y2) / n_y2
    result <- list(
        m_ub__y1 = m_ub__y1,
        pred_at_mesh__y1 = pred_at_mesh__y1,
        m_y2__y1 = m_y2__y1,
        Sigma_y2__y1 = Sigma_y2__y1,
        loo_error_list_squared = loo_error_list_squared,
        CRPS_list = CRPS_list,
        loo_error = loo_error,
        CRPS = CRPS,
        DSS = DSS,
        n_eff = n_eff,
        weights_normalized = weights_normalized,
        theta_sim_importance = theta_sim_importance,
        MAP = map_precip,
        MAP_par_non_log = convert_to_range_and_non_log(t(map_precip$par)),
        credible_interval = credible_interval
    )
    if (!is.null(path_data)) {
        saveRDS(result, file = path_data)
    }
    if (!is.null(path_prediction_txt)) {
        write.table(pred_at_mesh__y1, file = path_prediction_txt, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    return(result)
}

#' @title Leave out half testing
#' @description
#' Calculates the validation for the precipitation model using half the points for training and half for testing
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param location_indexes A vector of indices representing the training data. The testing data is the complement of the training data
#' @param data A list with x and y representing the observed data.
#' @param h_mesh A numeric vector representing the covariate at the mesh points
#' @param max_iterations_MAP Maximum number of iterations for MAP, by default 300
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (-0.5, c(0.1, 0.1), 0, -3, 0, 0)
#' @param n_weights The number of importance weights, by default 1000
#' @param path_data Path to save the result
#' @param path_prediction_txt Path to save the prediction
#' @return The calculated scores and estimates
#' @export

LOH_scores <- function(log_prior_theta, log_sigma_beta, mesh, location_indexes, data, h_mesh, max_iterations_MAP, theta0, n_weights, path_data = NULL) {
    observation_locations <- cbind(data$x[location_indexes], data$y[location_indexes])
    y1 <- data$precipitation[location_indexes] / 1000
    h1 <- data$altitude[location_indexes] / 1000
    A1 <- fm_basis(mesh, loc = observation_locations)
    y2 <- data$precipitation[-location_indexes] / 1000
    h2 <- data$altitude[-location_indexes] / 1000
    A2 <- fm_basis(mesh, loc = cbind(data$x[-location_indexes], data$y[-location_indexes]))

    LOH_score_2__1 <- LOM_scores(
        log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
        y1 = y1, A1 = A1, h1 = h1, y2 = y2, A2 = A2, h2 = h2, max_iterations_MAP = max_iterations_MAP,
        theta0 = theta0, n_weights = n_weights, h_mesh = h_mesh
    )
    LOH_score_1__2 <- LOM_scores(
        log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
        y1 = y2, A1 = A2, h1 = h2, y2 = y1, A2 = A1, h2 = h1, max_iterations_MAP = max_iterations_MAP,
        theta0 = theta0, n_weights = n_weights, h_mesh = h_mesh
    )
    # Takes average of each entry of scores and returns it
    loo_error <- (LOH_score_2__1$loo_error + LOH_score_1__2$loo_error) / 2
    CRPS <- (LOH_score_2__1$CRPS + LOH_score_1__2$CRPS) / 2
    DSS <- (LOH_score_2__1$DSS + LOH_score_1__2$DSS) / 2
    MAP_parameters_non_log <- (LOH_score_2__1$MAP_par_non_log + LOH_score_1__2$MAP_par_non_log) / 2
    credible_interval <- (LOH_score_2__1$credible_interval + LOH_score_1__2$credible_interval) / 2
    prediction_at_mesh <- (LOH_score_2__1$pred_at_mesh__y1 + LOH_score_1__2$pred_at_mesh__y1) / 2
    m_ub__y12 <- (LOH_score_2__1$m_ub__y1 + LOH_score_1__2$m_ub__y1) / 2
    result <- list(
        loo_error = loo_error,
        CRPS = CRPS,
        DSS = DSS,
        MAP_parameters_non_log = MAP_parameters_non_log,
        credible_interval = credible_interval,
        prediction_at_mesh = prediction_at_mesh,
        m_ub__y12 = m_ub__y12,
        LOH_score_2__1 = LOH_score_2__1,
        LOH_score_1__2 = LOH_score_1__2
    )
    if (!is.null(path_data)) {
        saveRDS(result, file = path_data)
    }
    return(result)
}

#' @title Leave out half loop
#' @description
#' Calculates the validation for the precipitation model using half the points for training and half for testing
#' @param log_prior_theta A function that calculates the log prior of theta = (log(kappa), v, log(sigma_u), log(sigma_epsilon))
#' @param log_sigma_beta The log variance of the linear effect beta
#' @param mesh The mesh
#' @param data A list with x and y representing the observed data
#' @param location_indexes_matrix A matrix. Each column are the indices of the data that is observed
#' @param max_iterations_MAP Maximum number of iterations for MAP, by default 300
#' @param theta0 Initial value for the parameters (log(kappa), v, log(sigma_u), log(sigma_epsilon)). By default, set to (-0.5, c(0.1, 0.1), 0, -3, 0, 0)
#' @param n_weights The number of importance weights, by default 1000
#' @param path_data Path to the data
#' @param path_prediction_txt Path to the prediction
#' @return The calculated average of the leave out half scores
#' @export

LOH_scores_loop <- function(log_prior_theta, log_sigma_beta, mesh, data, location_indexes_matrix, max_iterations_MAP = 300, theta0 = c(-0.5, c(0.1, 0.1), 0, -3, 0, 0), n_weights = 1000, h_mesh, path_data = NULL, path_prediction_txt = NULL) {
    n_loops <- dim(location_indexes_matrix)[2]
    observation_locations_list <- lapply(1:n_loops, function(i) {
        cbind(data$x[location_indexes_matrix[, i]], data$y[location_indexes_matrix[, i]])
    })

    A_list <- lapply(1:n_loops, function(i) {
        fm_basis(mesh, loc = observation_locations_list[[i]])
    })

    # Normalize precipitation and height to be approximately 1
    y_vector <- sapply(1:n_loops, function(i) {
        data$precipitation[location_indexes_matrix[, i]] / 1000
    })
    h_vector <- sapply(1:n_loops, function(i) {
        data$altitude[location_indexes_matrix[, i]] / 1000
    })
    n_loop <- dim(location_indexes_matrix)[2]
    n_y <- dim(location_indexes_matrix)[1]
    n_u <- nrow(mesh$loc)
    if (length(formals(log_prior_theta)) == 4) {
        d_theta <- 5
    } else if (length(formals(log_prior_theta)) == 3) {
        d_theta <- 3
    } else {
        stop("The prior function must have 3 or 4 arguments")
    }
    # Apply LOH_scores to each observation
    results <- lapply(1:n_loops, function(i) {
        LOH_scores(
            log_prior_theta = log_prior_theta, log_sigma_beta = log_sigma_beta, mesh = mesh,
            location_indexes = location_indexes_matrix[, i], data = data, h_mesh = h_mesh,
            max_iterations_MAP = max_iterations_MAP, theta0 = theta0, n_weights = n_weights
        )
    })
    # List of scores in each loop
    loo_error_list <- sapply(results, function(x) x$loo_error)
    CRPS_list <- sapply(results, function(x) x$CRPS)
    DSS_list <- sapply(results, function(x) x$DSS)
    # Get mean of scores
    loo_error <- mean(loo_error_list)
    CRPS <- mean(CRPS_list)
    DSS <- mean(DSS_list)
    # Get bootstrap of scores
    loo_error_bootstrap_CI <- bootstrap_q(loo_error_list, B = 1000)
    CRPS_bootstrap_CI <- bootstrap_q(CRPS_list, B = 1000)
    DSS_bootstrap_CI <- bootstrap_q(DSS_list, B = 1000)
    # Get mean map and credible interval
    MAP_non_log_par <- rowMeans(sapply(results, function(x) x$MAP))
    credible_interval <- matrix(rowMeans(sapply(results, function(x) x$credible_interval)), ncol = 2)
    # return the result
    result <- list(
        loo_error = loo_error,
        CRPS = CRPS,
        DSS = DSS,
        loo_error_list = loo_error_list,
        CRPS_list = CRPS_list,
        DSS_list = DSS_list,
        loo_error_bootstrap_CI = loo_error_bootstrap_CI,
        CRPS_bootstrap_CI = CRPS_bootstrap_CI,
        DSS_bootstrap_CI = DSS_bootstrap_CI,
        MAP_non_log_par = MAP_non_log_par,
        credible_interval = credible_interval
    )
    if (!is.null(path_data)) {
        saveRDS(result, file = path_data)
    }
    return(result)
}






#' @title DSS
#' @description
#' Calculates the Dawid-Sebastiani score
#' @param y The observed data
#' @param Q The posterior precision of the variable of interest
#' @param mu The posterior mean of the variable of interest
#' @return The Dawid-Sebastiani score (smaller is better)
#' @export

DSS <- function(mu, Q, y) {
    if (length(mu) == 1) {
        return(-log(Q) + t(y - mu) * Q * (y - mu))
    } else {
        return(-log(det(Q)) + t(y - mu) %*% Q %*% (y - mu))
    }
}


#' @title Bootstrap
#' @description
#' Given data s_1,...,s_n, calculates bootstrap estimate of the 95% confidence interval
#' @param data The data
#' @param B The number of bootstrap samples by default 10000
#' return The bootstrap estimate of the 95% confidence interval of the mean
#' @export

bootstrap_q <- function(data, B = 10000) {
    n <- length(data)
    boot_sample <- sample(data, size = n * B, replace = TRUE)
    boot_sample <- matrix(boot_sample, nrow = B)
    boot_means <- apply(boot_sample, 1, mean)
    boot_means <- sort(boot_means)
    lower <- boot_means[round(0.025 * B)]
    upper <- boot_means[round(0.975 * B)]
    data_mean <- mean(data)
    return(c(2 * data_mean - upper, 2 * data_mean - lower))
}

#' @title Bootstrap exchangeability
#' @description
#' Given two prediction scores {s_1,...,s_n}, {t_1,...,t_n}, calculates bootstrap estimate of the p-value for pairwise exchangeability
#' @param score1 The first score
#' @param score2 The second score
#' @param B The number of bootstrap samples by default 10000
#' return The bootstrap estimate of the p-value. The probability that model 1 is better than model 2
#' @export

bootstrap_pairwise_exchangeability <- function(data1, data2, B = 10000) {
    diff <- data1 - data2
    n <- length(diff)
    boot_sample <- sample(c(-1, 1), size = n * B, replace = TRUE) * sample(diff, size = n * B, replace = TRUE)
    boot_sample <- matrix(boot_sample, nrow = B)
    boot_means <- apply(boot_sample, 1, mean)
    data_mean <- mean(diff)
    boot_means <- sort(boot_means)
    p <- mean(boot_means > data_mean)
    return(p)
}
