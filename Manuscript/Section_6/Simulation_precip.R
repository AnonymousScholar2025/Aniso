# This file is used to obtain the results of the simulation
# comparing the various precipitation models
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
library(scoringRules)
library(progress)
document()

# # Get height first time
# observations_height <- getHeight(observation_locations)$hCov
# write.table(
#     observations_height,
#     "Precip/Matlab/sim_data/observations_height_sim.txt",
#     row.names = FALSE,
#     col.names = FALSE
# )


# Definition of priors ----------------------------------------------------
rho0 <- 10 # Controls the size of kappa in PC and non PC priors
a0 <- 10 # Controls the size of v in PC and non PC priors
sigma_u0 <- 3 # controls standard deviation of field
sigma_epsilon0 <- 3 # control standard deviation of noise
alpha <- 0.01 # quantile
log_sigma_beta <- log(10^2)

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
log_prior_theta_iso <- log_pc_prior_iso(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  rho0 = rho0,
  alpha = alpha
)


# Get locations and height of locations ---------------------------------------------------------------
load("Manuscript/Section_6/Precip/data/data.RData")
mesh <- mesh
h_mesh <- scan("Manuscript/Section_6/Precip/data/mesh_elevation.txt")
border <- read.table("Manuscript/Section_6/Precip/Matlab/delBorder.txt")
observation_locations <- as.matrix(read.csv("Manuscript/Section_6/Precip/Matlab/sim_data/observation_locations_sim.txt"))
observations_height <- scan("Manuscript/Section_6/Precip/Matlab/sim_data/observations_height_sim.txt")
n_loc_max <- length(observations_height) # Max number of locations

# Get simulated data ------------------------------------------------------
simulation_type <- "ANISO" # "ANISO" or "ISO"

if (simulation_type == "ANISO") {
  LOO <- readRDS("Manuscript/Section_6/Precip/Results/LOO_pc_20000.RData")
  theta0 <- LOO$MAP$par
  v <- theta0[2:3]
  log_sigma_u <- theta0[4]
  sigma_epsilon <- exp(theta0[5])
} else if (simulation_type == "ISO") {
  LOO <- readRDS("Manuscript/Section_6/Precip/Results/LOO_iso_10000.RData")
  theta0 <- LOO$MAP$par
  theta0 <- c(theta0[1], 10^-6, 0, theta0[2], theta0[3])
  v <- c(0, 0)
  log_sigma_u <- theta0[4]
  sigma_epsilon <- exp(theta0[5])
}
beta0 <- LOO$m_ub__y[1]
beta1 <- LOO$m_ub__y[2]
theta0_non_log <- convert_to_range_and_non_log(t(theta0))
kappa <- exp(theta0[1])
aniso <- list(rep(kappa, mesh$n), matrix(v, mesh$n, 2, byrow = TRUE))

u <- fm_aniso_basis_weights_sample(
  x = mesh,
  aniso = aniso,
  log_sigma = log_sigma_u
)

A <- fmesher::fm_basis(mesh, loc = observation_locations)
noise <- rnorm(n_loc_max, 0, sigma_epsilon)
precip <- 1000 * (A %*% u + noise +
  beta0) + beta1 * observations_height
data_sim <- data.frame(
  x = observation_locations[, 1],
  y = observation_locations[, 2],
  precipitation = as.vector(precip),
  altitude = observations_height
)

#Save data
save(data_sim,
  file = paste0(
    "Manuscript/Section_6/Precip/Matlab/sim_data/data_sim_",
    simulation_type,
    ".RData"
  )
)
# Load data
load(paste0(
  "Manuscript/Section_6/Precip/Matlab/sim_data/data_sim_",
  simulation_type,
  ".RData"
))


# Plotting data -----------------------------------------------------------

# Create dataframes for ggplot
field_mesh <- data.frame(
  x = mesh$loc[, 1],
  y = mesh$loc[, 2],
  u = as.vector(u)
)
field_on_obs_locations <- data.frame(
  x = observation_locations[, 1],
  y = observation_locations[, 2],
  u = as.vector(A %*% u)
)

# Plot field on mesh
# ggplot() +
#   geom_path(aes(x = border$V1, y = border$V2)) +
#   geom_point(data = field_mesh, aes(x = x, y = y, color = u), size = 1) +
#   scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
#   theme_void()

# Plot field on observation locations
ggplot() +
  geom_path(aes(x = border$V1, y = border$V2)) +
  geom_point(data = field_on_obs_locations, aes(x = x, y = y, color = u), size = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  # ggtitle("Simulated precipitation") +
  xlab("Easting (km)") +
  ylab("Northing (km)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 25)) +
  xlim(-170, 630) +
  ylim(6350, 7300) +
  guides(color = guide_colorbar(title = NULL, barheight = 30, barwidth = 1))
ggsave(
  paste0(
    "Manuscript/Section_6/Precip/figures/simulated_field_",
    simulation_type,
    ".pdf"
  ),
  width = 10,
  height = 10
)

ggplot() +
  geom_path(aes(x = border$V1, y = border$V2)) +
  geom_point(data = data_sim, aes(x = x, y = y, color = precipitation / 1000), size = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  # ggtitle("Simulated precipitation") +
  xlab("Easting (km)") +
  ylab("Northing (km)") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 25)) +
  xlim(-170, 630) +
  ylim(6350, 7300) +
  guides(color = guide_colorbar(title = NULL, barheight = 30, barwidth = 1))
ggsave(
  paste0(
    "Manuscript/Section_6/Precip/figures/simulated_data_",
    simulation_type,
    ".pdf"
  ),
  width = 10,
  height = 10
)
# Plot height
ggplot() +
  geom_path(aes(x = border$V1, y = border$V2)) +
  geom_point(data = data_sim, aes(x = x, y = y, color = altitude / 1000), size = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_void()


# Parameters of simulation ------------------------------------------------
n_locations <- 25 # Number of locations to observe (300,400,...,1000)




# Unchanging parameters
n_loops <- 100 # Number of times to run the simulation (unchanging)
n_weights <- 100 # Number of weights to use in importance sampling (unchanging)

max_iterations_MAP <- 600 # max iterations for the optimization
# Locations observed for each loop
set.seed(123)
location_indexes_matrix <- replicate(n_loops, sample.int(length(data_sim$x))[1:n_locations])



# Calculation of LOO ------------------------------------------------------

LOO_pc <- LOO_CV_precipitation_loop(
  log_prior_theta = log_pc_prior,
  log_sigma_beta = log_sigma_beta,
  mesh = mesh,
  h_mesh = h_mesh,
  data = data_sim,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Manuscript/Section_6/Precip/data_sim/",
    simulation_type,
    "/LOO_pc_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh_sim/",
    simulation_type,
    "/prediction_on_mesh_pc_w=",
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
  data = data_sim,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Manuscript/Section_6/Precip/data_sim/",
    simulation_type,
    "/LOO_EG_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh_sim/",
    simulation_type,
    "/prediction_on_mesh_EG_w=",
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
  data = data_sim,
  location_indexes_matrix = location_indexes_matrix,
  max_iterations_MAP = max_iterations_MAP,
  theta0 = theta0,
  n_weights = n_weights,
  path_data = paste0(
    "Precip/data_sim/",
    simulation_type,
    "/LOO_iso_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".RData"
  ),
  path_prediction_txt = paste0(
    "Manuscript/Section_6/Precip/Matlab/Prediction_on_mesh_sim/",
    simulation_type,
    "/prediction_on_mesh_iso_w=",
    n_weights,
    "_loc=",
    n_locations,
    "_n_loops=",
    n_loops,
    ".txt"
  )
)


# Extracting interesting data from LOO ---------------------------------------------

# n_locations <- 2000
# simulation_type <- "ISO"
# n_weights <- 100
# n_loops <- n_loops
# # Results
# LOO_pc <- readRDS(
#   paste0(
#     "Manuscript/Section_6/Precip/data_sim/",
#     simulation_type,
#     "/LOO_pc_w=",
#     n_weights,
#     "_loc=",
#     n_locations,
#     "_n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# LOO_pc$interval

# LOO_EG <- readRDS(
#   paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOO_EG_w=",
#     n_weights,
#     "_loc=",
#     n_locations,
#     "_n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# LOO_iso <- readRDS(
#   paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOO_iso_w=",
#     n_weights,
#     "_loc=",
#     n_locations,
#     "_n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# # MAPS
# bootstrap_pairwise_exchangeability(LOO_iso$loo_error_list, LOO_pc$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$loo_error_list, LOO_pc$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_iso$DSS_list, LOO_pc$DSS_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$DSS_list, LOO_pc$DSS_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_iso$CRPS_list, LOO_pc$CRPS_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$CRPS_list, LOO_pc$CRPS_list, 10000)
# MAPS_non_log <- list(
#   pc = LOO_pc$MAP_non_log_par,
#   EG = LOO_EG$MAP_non_log_par,
#   iso = LOO_iso$MAP_non_log_par
# )

# loo_error <- list(
#   pc = LOO_pc$loo_error,
#   EG = LOO_EG$loo_error,
#   iso = LOO_iso$loo_error
# )

# CRPS <- list(
#   pc = LOO_pc$CRPS,
#   EG = LOO_EG$CRPS,
#   iso = LOO_iso$CRPS
# )

# DSS <- list(
#   pc = LOO_pc$DSS,
#   EG = LOO_EG$DSS,
#   iso = LOO_iso$DSS
# )

# loo_error_bootstrap_CI <- list(
#   pc = LOO_pc$loo_error_bootstrap_CI,
#   EG = LOO_EG$loo_error_bootstrap_CI,
#   iso = LOO_iso$loo_error_bootstrap_CI
# )
# bootstrap_pairwise_exchangeability(LOO_iso$loo_error_list, LOO_pc$loo_error_list, 10000)
# bootstrap_pairwise_exchangeability(LOO_EG$loo_error_list, LOO_pc$loo_error_list, 10000)




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
# #            path_data = paste0("Precip/data_sim/LOM_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"),
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
# #                     path_data = paste0("Precip/data_sim/LOM_pc_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".RData"),
# #                     path_prediction_txt = paste0("Precip/Matlab/Prediction_on_mesh/prediction_on_mesh_iso_w=", n_weights, "_loc=", n_locations, "_n_loops=", n_loops, ".txt")
# # )




# Leave out half loop -----------------------------------------------------


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
#   path_data = paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_pc_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
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
#   path_data = paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_EG_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
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
#   path_data = paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_iso_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
#
# # Extracting scores from LOH loop
#
# n_weights <- n_weights
# n_loops <- n_loops
#
# LOH_pc_loop <- readRDS(
#   paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_pc_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# LOH_EG_loop <- readRDS(
#   paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_EG_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# LOH_iso_loop <- readRDS(
#   paste0(
#     "Precip/data_sim/",
#     simulation_type,
#     "/LOH_iso_w=",
#     n_weights,
#     "n_loops=",
#     n_loops,
#     ".RData"
#   )
# )
# bootstrap_pairwise_exchangeability(LOH_iso_loop$loo_error_list,
#                                    LOH_pc_loop$loo_error_list,
#                                    10000)
# bootstrap_pairwise_exchangeability(LOH_EG_loop$loo_error_list,
#                                    LOH_pc_loop$loo_error_list,
#                                    10000)
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
# CRPS_half <- list(pc = LOH_pc_loop$CRPS,
#                   EG = LOH_EG_loop$CRPS,
#                   iso = LOH_iso_loop$CRPS)
#
#
# DSS_half <- list(pc = LOH_pc_loop$DSS,
#                  EG = LOH_EG_loop$DSS,
#                  iso = LOH_iso_loop$DSS)
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
