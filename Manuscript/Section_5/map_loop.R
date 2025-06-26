library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(future)
library(future.apply)
library(dplyr)
library(tidyr)
library(loo)
library(progress)
document()

# Defining the random seed
set.seed(123)

# Distributions to sample true parameters from ----------------------------------------------------------
prior_distributions <- c("pc", "EG", "uniform", "beta")
true_parameter_distribution <- prior_distributions[4]

# Hyerparameters ----------------------------------------------------------
rho0 <- 1 # Controls the size of kappa
a0 <- 2 # Controls the size of v
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 2 # control standard deviation of noise
sigma0 <- 1.5 # Controls the size of v in non PC priors

L <- 10

shape <- 1.1
# beta and uniform shape is smaller if simulated from to not get extreme values
if (true_parameter_distribution == "beta") {
  width_beta <- 2
} else {
  width_beta <- 20
}
if (true_parameter_distribution == "uniform") {
  width_uniform <- 2
} else {
  width_uniform <- Inf
}

alpha <- 0.01 # Defines the quantile

m_u <- 0 # Setting mean of the field



# Priors -------------------------------------------------------------------

log_pc_prior <- log_pc_prior_quantile(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  alpha = alpha
)

log_EG_prior <- log_EG_prior_quantile(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  alpha = alpha
)



log_beta_prior <- log_prior_beta(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  L = L,
  shape = shape,
  width_support_factor = width_beta
)
log_uniform_prior <- log_prior_uniform(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  L = L,
  width_support_factor = Inf
)
log_priors <- list(
  pc = log_pc_prior,
  EG = log_EG_prior,
  uniform = log_uniform_prior,
  beta = log_beta_prior
)


prior_types <- setNames(as.list(names(log_priors)), names(log_priors))
approximation_types <- list("Gaussian_median", "importance", "importance_smoothed")
approximation_types <- setNames(approximation_types, approximation_types)


# Mesh and observation --------------------------------------------------------------------
library(sf)
boundary_sf <- st_sfc(st_polygon(list(rbind(
  c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)
))))
boundary <- fmesher::fm_as_segm(boundary_sf)
mesh <- fmesher::fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 3))
nodes <- mesh$loc
n <- mesh$n
# plot(mesh)
n_observations <- 15
observations <- L * matrix(runif(n_observations * 2), ncol = 2) # Points where the field is observed
A <- fmesher::fm_basis(mesh, loc = observations)

number_of_loops <- 600 # number of iterations
number_of_weights <- 1000
maxit_MAP <- 600
credible_level <- 0.05
results <- vector("list", number_of_loops) # Pre-allocates a list for number_of_loops iterations
set.seed(123)
pb <- progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = number_of_loops,
  complete = "=",
  incomplete = "-",
  current = ">",
  clear = TRUE,
  width = 100
)
for (i in 1:number_of_loops) {
  # start_time <- Sys.time()
  tryCatch(
    {
      # Simulate parameters from PC prior
      if (true_parameter_distribution == "pc") {
        true_params <- sim_theta_pc_quantile(
          alpha = alpha,
          sigma_u0 = sigma_u0,
          sigma_epsilon0 = sigma_epsilon0,
          a0 = a0,
          rho0 = rho0,
          m = 1
        )
      } else if (true_parameter_distribution == "EG") {
        true_params <- sim_not_pc(
          alpha = alpha,
          sigma_u0 = sigma_u0,
          sigma_epsilon0 = sigma_epsilon0,
          a0 = a0,
          rho0 = rho0,
          m = 1
        )
      } else if (true_parameter_distribution == "uniform") {
        true_params <- sim_theta_uniform(
          alpha = alpha,
          sigma_u0 = sigma_u0,
          sigma_epsilon0 = sigma_epsilon0,
          a0 = a0,
          rho0 = rho0,
          L = L,
          width_support_factor = width_uniform
        )
      } else if (true_parameter_distribution == "beta") {
        true_params <- sim_theta_beta(
          alpha = alpha,
          sigma_u0 = sigma_u0,
          sigma_epsilon0 = sigma_epsilon0,
          a0 = a0,
          rho0 = rho0,
          L = L,
          shape = shape,
          width_support_factor = width_beta
        )
      }
      log_kappa <- true_params$log_kappa
      kappa <- exp(log_kappa)
      v <- true_params$v
      log_sigma_u <- true_params$log_sigma_u
      log_sigma_epsilon <- true_params$log_sigma_epsilon
      # aniso <- list(rep(kappa, n), matrix(v, n, 2)) old version was incorrect
      aniso <- list(rep(kappa, n), matrix(v, mesh$n, 2, byrow = TRUE))

      # Sample from noisy data
      x <- fm_aniso_basis_weights_sample(
        x = mesh,
        aniso = aniso,
        log_sigma = log_sigma_u
      )
      y <- A %*% x + rnorm(n_observations, 0, exp(log_sigma_epsilon))

      maps <- lapply(log_priors, function(log_prior) {
        MAP_prior(
          log_prior = log_prior,
          mesh = mesh,
          y = y,
          A = A,
          m_u = m_u,
          max_iterations = maxit_MAP,
          theta0 = unlist(true_params)
        )
      })

      # Gaussian_median approximations
      mus_Gaussian_median <- lapply(maps, function(map) {
        map$par
      })
      Qs_Gaussian_median <- lapply(maps, function(map) {
        -map$hessian
      })
      Covariances_Gaussian_median <- lapply(Qs_Gaussian_median, function(Q) {
        solve(Q)
      })
      std_dev_estimates_Gaussian_median <- lapply(Covariances_Gaussian_median, function(Covariance) {
        sqrt(diag(Covariance))
      })

      m_u <- 0

      log_posteriors <- lapply(log_priors, function(log_prior) {
        function(theta) {
          log_kappa <- theta[1]
          v <- theta[2:3]
          log_sigma_u <- theta[4]
          log_sigma_epsilon <- theta[5]
          log_posterior_prior(
            log_prior = log_prior,
            mesh = mesh,
            log_kappa = log_kappa,
            v = v,
            log_sigma_epsilon = log_sigma_epsilon,
            log_sigma_u = log_sigma_u,
            y = y,
            A = A,
            m_u = m_u
          )
        }
      })

      # Importance sampling
      importances <- lapply(prior_types, function(prior_type) {
        log_posterior <- log_posteriors[[prior_type]]
        mu_Gaussian_median <- mus_Gaussian_median[[prior_type]]
        Q_Gaussian_median <- Qs_Gaussian_median[[prior_type]]
        log_unnormalized_importance_weights_and_integrals(
          log_posterior_density = log_posterior,
          mu_Gaussian_median = mu_Gaussian_median,
          Q_Gaussian_median = Q_Gaussian_median,
          n_weights = number_of_weights,
          q = credible_level,
          true_params = unlist(true_params)
        )
      })

      # CIs
      credible_intervals <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          importances[[prior_type]][[paste0("credible_intervals_", approximation_type)]]
        })
      })

      true_parameter_is_within_CI <- lapply(prior_types, function(prior_type) {
        lapply(approximation_types, function(approximation_type) {
          parameter_within_credible_intervals(true_params, credible_intervals[[prior_type]][[approximation_type]])
        })
      })


      # Accumulate results
      results_accumulator <- function(prior_type) {
        # Return a list of the calculated values
        list(
          MAP_estimate = maps[[prior_type]]$par,
          MAP_value = maps[[prior_type]]$value,
          convergence = maps[[prior_type]]$convergence,
          distance_vector = abs(maps[[prior_type]]$par - unlist(true_params)),
          covariance_estimate = Covariances_Gaussian_median[[prior_type]],
          std_dev_estimates_Gaussian_median = std_dev_estimates_Gaussian_median[[prior_type]],
          credible_intervals = lapply(approximation_types, function(approximation_type) {
            credible_intervals[[prior_type]][[approximation_type]]
          }),
          true_parameter_within_c_interval = lapply(approximation_types, function(approximation_type) {
            true_parameter_is_within_CI[[prior_type]][[approximation_type]]
          }),
          importance = importances[[prior_type]]
        )
      }


      # Store results
      partial_results <- lapply(prior_types, results_accumulator)
      results[[i]] <- c(
        list(true_params = true_params),
        partial_results,
        prior_types
      )
    },
    error = function(e) {
      e
    }
  )
  pb$tick()
}
# Eliminates NULL results
not_null_indices <- sapply(results, function(x) {
  !is.null(x$pc$importance$n_eff)
})
results <- results[not_null_indices]
# Results obtained simulating parameters from PC priors and using a mesh size of 1, 15 observations, 200 iterations, 5000 weights, a credible level of 0.05 a width of uniform =inf and for beta a multiplier of 20.
saveRDS(results, paste0(
  "Manuscript/Section_5/Simulation_results/results_", true_parameter_distribution,
  "_n_obs=",
  n_observations,
  "_n_loops=",
  number_of_loops,
  "_n_weights=",
  number_of_weights,
  "_wu=",
  width_uniform,
  "_wb=",
  width_beta,
  ".rds"
))
