# Visualizing posteriors and priors
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
load("Precipitation/data.RData")
mesh <- mesh
h_mesh <- scan("Precip/data/mesh_elevation.txt")
border <- read.table("Precip/Matlab/delBorder.txt")
observation_locations <- as.matrix(read.csv("Precip/Matlab/sim_data/observation_locations_sim.txt"))
observations_height <- scan("Precip/Matlab/sim_data/observations_height_sim.txt")
n_loc_max <- length(observations_height) # Max number of locations

# Get simulated data ------------------------------------------------------

theta0 <- c(-4.26, 3, 0, -0.45, 2 * -1.953)
theta0_non_log <- convert_to_range_and_non_log(t(theta0))
kappa <- exp(theta0[1])
v <- theta0[2:3]
log_sigma_u <- theta0[4]
sigma_epsilon <- exp(theta0[5])
# aniso <- list(rep(kappa, mesh$n), matrix(v, mesh$n, 2, byrow = TRUE))
aniso <- list(rep(kappa, mesh$n), matrix(v, mesh$n, 2))

u_bad <- fm_aniso_basis_weights_sample(
  x = mesh,
  aniso = list(rep(kappa, mesh$n), matrix(v, mesh$n, 2)),
  log_sigma = log_sigma_u
)
u_good <- fm_aniso_basis_weights_sample(
  x = mesh,
  aniso = list(rep(kappa, mesh$n), matrix(v, mesh$n, 2, byrow = TRUE)),
  log_sigma = log_sigma_u
)

A <- fmesher::fm_basis(mesh, loc = observation_locations)
noise <- rnorm(n_loc_max, 0, sigma_epsilon)
y_good <- A %*% u_good + noise
y_bad <- A %*% u_bad + noise

# Plotting data -----------------------------------------------------------

# Create dataframes for ggplot
field_on_obs_locations <- data.frame(
  x = observation_locations[, 1],
  y = observation_locations[, 2],
  u_good = as.vector(A %*% u_good + noise),
  u_bad = as.vector(A %*% u_bad + noise)
)

# Plot field on observation locations
ggplot() +
  geom_path(aes(x = border$V1, y = border$V2)) +
  geom_point(data = field_on_obs_locations, aes(x = x, y = y, color = u_good), size = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_void()

ggplot() +
  geom_path(aes(x = border$V1, y = border$V2)) +
  geom_point(data = field_on_obs_locations, aes(x = x, y = y, color = u_bad), size = 1) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_void()

n_locations <- 1000 # Number of locations to observe (300,400,...,1000)

mp_good <- MAP_prior(
  log_prior = log_pc_prior,
  mesh = mesh,
  y = y_good,
  A = A,
  m_u = 0,
  max_iterations = 600,
  theta0 = theta0
)
mp_bad <- MAP_prior(
  log_prior = log_pc_prior,
  mesh = mesh,
  y = y_bad,
  A = A,
  m_u = 0,
  max_iterations = 600,
  theta0 = theta0
)


log_post <- function(theta) {
  log_kappa <- theta[1]
  v <- theta[2:3]
  log_sigma_u <- theta[4]
  log_sigma_epsilon <- theta[5]
  log_posterior_prior(
    log_prior = log_pc_prior,
    mesh = mesh, log_kappa = log_kappa, v = v,
    log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon,
    y = y, A = A, m_u = 0
  )
}
mp_good$par
mp_bad$par
theta0
log_post(theta0)
log_post(c(-4.26, 3.00, 3.00, -0.45, -4))
log_post(mp$par)
