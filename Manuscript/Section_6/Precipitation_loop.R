#This file is used to compare the various precipitation models under varying
#ammounts of observed data
library(devtools)
library(SPDEaniso)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(loo)
library(foreach)
library(doParallel)
install.packages("scoringRules")
library(scoringRules)
library(progress)
#devtools::document(pkg = "package")

# Defines the upper bounds for the quantiles
rho0 <- 10 # Controls the size of kappa in PC and non PC priors
a0 <- 10 # Controls the size of v in PC and non PC priors
sigma_u0 <- 3 # controls standard deviation of field
sigma_epsilon0 <- 3 # control standard deviation of noise
alpha <- 0.01 # quantile
log_sigma_beta <- log(10^2)

log_pc_prior <- log_pc_prior_quantile(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)

log_EG_prior <- log_EG_prior_quantile(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  a0 = a0, rho0 = rho0, alpha = alpha
)
log_prior_theta_iso <- log_pc_prior_iso(
  sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
  rho0 = rho0, alpha = alpha
)


# Load and extract data
load("Manuscript/Section_6/Precip/data/data.RData")
mesh <- fm_as_mesh_2d(mesh)
h_mesh <- scan("Manuscript/Section_6/Precip/data/mesh_elevation.txt")

# Parameters of loop
n_locations <- 100 # Number of observations to observe
n_loops <- 100 # Number of times to run the simulation
max_iterations_MAP <- 600 # max iterations for the optimization
n_weights <- 100 # Number of weights to use in importance sampling

# Locations observed for each loop
set.seed(123)
location_indexes_matrix <- replicate(n_loops, sample.int(length(data$x))[1:n_locations])

# Simulates the starting point for optimization from the pc_prior.
theta0 <- unlist(sim_theta_pc_quantile(
  alpha = alpha, sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1
))


# Calculation of LOO ------------------------------------------------------

LOO_pc <- LOO_CV_precipitation_loop(
  log_prior_theta = log_pc_prior,
  log_sigma_beta = log_sigma_beta,
  mesh = mesh,
  h_mesh = h_mesh,
  data = data,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Manuscript/Section_6/Precip/Results/LOO_pc_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_pc_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".txt"
  )
)
LOO_EG <- LOO_CV_precipitation_loop(
  log_prior_theta = log_EG_prior,
  log_sigma_beta = log_sigma_beta,
  mesh = mesh,
  h_mesh = h_mesh,
  data = data,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Manuscript/Section_6/Precip/Results/LOO_EG_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_EG_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".txt"
  )
)
LOO_iso <- LOO_CV_precipitation_loop(
  log_prior_theta = log_prior_theta_iso,
  log_sigma_beta = log_sigma_beta,
  mesh = mesh,
  h_mesh = h_mesh,
  data = data,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Manuscript/Section_6/Precip/Results/LOO_iso_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_iso_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".txt"
  )
)


# # Extracting interesting data from LOO ---------------------------------------------

# # n_locations <- n_locations
# n_weights <- n_weights
# n_loops <- n_loops
# # Results
# LOO_pc <- readRDS(paste0("Manuscript/Section_6/Precip/Results/LOO_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"))
# LOO_EG <- readRDS(paste0("Manuscript/Section_6/Precip/Results/LOO_EG_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"))
# LOO_iso <- readRDS(paste0("Manuscript/Section_6/Precip/Results/LOO_iso_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"))
# # MAPS
# bootstrap_pairwise_exchangeability(LOO_iso$loo_error_list, LOO_pc$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$loo_error_list, LOO_pc$loo_error_list, 10000)
#
# MAPS_non_log <- list(
#   pc = LOO_pc$MAP_non_log_par,
#   EG = LOO_EG$MAP_non_log_par,
#   iso = LOO_iso$MAP_non_log_par
# )
#
# loo_error <- list(
#   pc = LOO_pc$loo_error,
#   EG = LOO_EG$loo_error,
#   iso = LOO_iso$loo_error
# )
#
# CRPS <- list(
#   pc = LOO_pc$CRPS,
#   EG = LOO_EG$CRPS,
#   iso = LOO_iso$CRPS
# )
#
# DSS <- list(
#   pc = LOO_pc$DSS,
#   EG = LOO_EG$DSS,
#   iso = LOO_iso$DSS
# )
#
# loo_error_bootstrap_CI <- list(
#   pc = LOO_pc$loo_error_bootstrap_CI,
#   EG = LOO_EG$loo_error_bootstrap_CI,
#   iso = LOO_iso$loo_error_bootstrap_CI
# )
# bootstrap_pairwise_exchangeability(LOO_iso$loo_error_list, LOO_pc$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$loo_error_list, LOO_pc$loo_error_list, 10000)
#
#
#
# # Leave out half_loop
# set.seed(123)
# n_half <- round(length(data$precipitation) / 2)
# n_loops <- n_loops
# n_weights <- n_weights
# location_indexes_matrix <- replicate(n_loops, sample.int(length(data$x))[1:n_half])
# LOH_pc_loop <- LOH_scores_loop(
#   log_prior_theta = log_pc_prior,
#   log_sigma_beta = log_sigma_beta,
#   mesh = mesh,
#   location_indexes_matrix = location_indexes_matrix,
#   data = data,
#   h_mesh = h_mesh,
#   max_iterations_MAP = max_iterations_MAP,
#   theta0 = theta0,
#   n_weights = n_weights,
#   path_data = paste0("Precip/data/LOH_pc_w=", n_weights, "n_loops=", n_loops, ".RData")
# )
#
# LOH_EG_loop <- LOH_scores_loop(
#   log_prior_theta = log_EG_prior,
#   log_sigma_beta = log_sigma_beta,
#   mesh = mesh,
#   location_indexes_matrix = location_indexes_matrix,
#   data = data,
#   h_mesh = h_mesh,
#   max_iterations_MAP = max_iterations_MAP,
#   theta0 = theta0,
#   n_weights = n_weights,
#   path_data = paste0("Precip/data/LOH_EG_w=", n_weights, "n_loops=", n_loops, ".RData")
# )
#
# LOH_iso_loop <- LOH_scores_loop(
#   log_prior_theta = log_prior_theta_iso,
#   log_sigma_beta = log_sigma_beta,
#   mesh = mesh,
#   location_indexes_matrix = location_indexes_matrix,
#   data = data,
#   h_mesh = h_mesh,
#   max_iterations_MAP = max_iterations_MAP,
#   theta0 = theta0,
#   n_weights = n_weights,
#   path_data = paste0("Precip/data/LOH_iso_w=", n_weights, "n_loops=", n_loops, ".RData")
# )
#
# # Extracting scores from LOH loop
#
# n_weights <- n_weights
# n_loops <- n_loops
#
# LOH_pc_loop <- readRDS(paste0("Precip/data/LOH_pc_w=", n_weights, "n_loops=", n_loops, ".RData"))
# LOH_EG_loop <- readRDS(paste0("Precip/data/LOH_EG_w=", n_weights, "n_loops=", n_loops, ".RData"))
# LOH_iso_loop <- readRDS(paste0("Precip/data/LOH_iso_w=", n_weights, "n_loops=", n_loops, ".RData"))
#
# bootstrap_pairwise_exchangeability(LOH_iso_loop$loo_error_list, LOH_pc_loop$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_EG_loop$loo_error_list, LOH_pc_loop$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_iso_loop$CRPS_list, LOH_pc_loop$CRPS_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_EG_loop$CRPS_list, LOH_pc_loop$CRPS_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_iso_loop$DSS_list, LOH_pc_loop$DSS_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_iso_loop$DSS_list, LOH_EG_loop$DSS_list, 10000)
# bootstrap_pairwise_exchangeability(LOH_EG_loop$DSS_list, LOH_pc_loop$DSS_list, 10000)
#
#
# MAPS_non_log_half <- list(
#   pc = LOH_pc_loop$MAP_non_log_par,
#   EG = LOH_EG_loop$MAP_non_log_par,
#   iso = LOH_iso_loop$MAP_non_log_par
# )
#
# loo_error_half <- list(
#   pc = LOH_pc_loop$loo_error,
#   EG = LOH_EG_loop$loo_error,
#   iso = LOH_iso_loop$loo_error
# )
#
# CRPS_half <- list(
#   pc = LOH_pc_loop$CRPS,
#   EG = LOH_EG_loop$CRPS,
#   iso = LOH_iso_loop$CRPS
# )
#
#
# DSS_half <- list(
#   pc = LOH_pc_loop$DSS,
#   EG = LOH_EG_loop$DSS,
#   iso = LOH_iso_loop$DSS
# )
#
# loo_error_half_bootstrap_CI <- list(
#   pc = LOH_pc_loop$loo_error_bootstrap_CI,
#   EG = LOH_EG_loop$loo_error_bootstrap_CI,
#   iso = LOH_iso_loop$loo_error_bootstrap_CI
# )
#
# CRPS_half_bootstrap_CI <- list(
#   pc = LOH_pc_loop$CRPS_bootstrap_CI,
#   EG = LOH_EG_loop$CRPS_bootstrap_CI,
#   iso = LOH_iso_loop$CRPS_bootstrap_CI
# )
#
# DSS_half_bootstrap_CI <- list(
#   pc = LOH_pc_loop$DSS_bootstrap_CI,
#   EG = LOH_EG_loop$DSS_bootstrap_CI,
#   iso = LOH_iso_loop$DSS_bootstrap_CI
# )




#
# # Leave out half ----------------------------------------------------------
# set.seed(321)
# n_half <- round(length(data$precipitation)/2)
# n_weights <- 10
# location_indexes <- sample.int(length(data$x))[1:n_half]
# # LOH_pc<- LOH_scores(log_prior_theta = log_pc_prior,
# #            log_sigma_beta = log_sigma_beta,
# #            mesh = mesh,
# #            location_indexes = location_indexes,
# #            data = data,
# #            h_mesh = h_mesh,
# #            max_iterations_MAP = max_iterations_MAP,
# #            theta0 = theta0,
# #            n_weights = n_weights,
# #            path_data = paste0("Precip/data/LOH_pc_w=", n_weights, ".RData")
# # )
# # LOH_EG<- LOH_scores(log_prior_theta = log_EG_prior,
# #                     log_sigma_beta = log_sigma_beta,
# #                     mesh = mesh,
# #                     location_indexes = location_indexes,
# #                     data = data,
# #                     h_mesh = h_mesh,
# #                     max_iterations_MAP = max_iterations_MAP,
# #                     theta0 = theta0,
# #                     n_weights = n_weights,
# #                     path_data = paste0("Precip/data/LOH_EG_w=", n_weights, ".RData")
# # )
# LOH_iso<- LOH_scores(log_prior_theta = log_prior_theta_iso,
#                     log_sigma_beta = log_sigma_beta,
#                     mesh = mesh,
#                     location_indexes = location_indexes,
#                     data = data,
#                     h_mesh = h_mesh,
#                     max_iterations_MAP = max_iterations_MAP,
#                     theta0 = theta0,
#                     n_weights = n_weights,
#                     path_data = paste0("Precip/data/LOH_iso_w=", n_weights, ".RData")
# )
#
# #Extracting scores from LOH
#
# MAPS_non_log_half <- list(
#   pc = LOH_pc$MAP_parameters_non_log,
#   EG = LOH_EG$MAP_parameters_non_log,
#   iso = LOH_iso$MAP_parameters_non_log)
#
# loo_error_half <- list(
#   pc = LOH_pc$loo_error,
#   EG = LOH_EG$loo_error,
#   iso = LOH_iso$loo_error
# )
#
# crps_half <- list(
#   pc = LOH_pc$CRPS,
#   EG = LOH_EG$CRPS,
#   iso = LOH_iso$CRPS
# )
#
# dss_half <- list(
#   pc = LOH_pc$DSS,
#   EG = LOH_EG$DSS,
#   iso = LOH_iso$DSS
# )
# LOH_iso$loo_error
# LOH_iso$CRPS
# LOH_iso$DSS
# # Test LOM
# # set.seed(123)
# # n_locations <- length(data$precipitation) -150
# # n_weights <- 2
# # location_indexes <- sample.int(length(data$x))[1:n_locations]
# # y1 <- data$precipitation[location_indexes] / 1000
# # A1 <- fm_basis(mesh, loc = cbind(data$x[location_indexes], data$y[location_indexes]))
# # h1 <- data$altitude[location_indexes] / 1000
# # y2 <- data$precipitation[-location_indexes] / 1000
# # A2 <- fm_basis(mesh, loc = cbind(data$x[-location_indexes], data$y[-location_indexes]))
# # h2 <- data$altitude[-location_indexes] / 1000
# #
# #
# # LOM_pc<- LOM_scores(log_prior_theta = log_pc_prior,
# #            log_sigma_beta = log_sigma_beta,
# #            mesh = mesh,
# #            y1 = y1,
# #            A1 = A1,
# #            h1 = h1,
# #            y2 = y2,
# #            A2 = A2,
# #            h2 = h2,
# #            max_iterations_MAP = max_iterations_MAP,
# #            theta0 = theta0,
# #            n_weights = n_weights,
# #            h_mesh = h_mesh,
# #            path_data = paste0("Precip/data/LOM_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"),
# #            path_prediction_txt = paste0("Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".txt")
# # )
# #
# # LOM_iso<- LOM_scores(log_prior_theta = log_prior_theta_iso,
# #                     log_sigma_beta = log_sigma_beta,
# #                     mesh = mesh,
# #                     y1 = y1,
# #                     A1 = A1,
# #                     h1 = h1,
# #                     y2 = y2,
# #                     A2 = A2,
# #                     h2 = h2,
# #                     max_iterations_MAP = max_iterations_MAP,
# #                     theta0 = theta0,
# #                     n_weights = n_weights,
# #                     h_mesh = h_mesh,
# #                     path_data = paste0("Precip/data/LOM_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"),
# #                     path_prediction_txt = paste0("Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_iso_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".txt")
# # )
