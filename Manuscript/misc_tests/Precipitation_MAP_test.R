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
library(sf)
# document()#

load("Precipitation/data.RData")
# locations are data$x, data$y we turn into a matrix
mesh <- mesh # Mesh of the domain, is inla.mesh not fm_mesh_2d, so unsure how to plot
mesh_points <- as.data.frame(mesh$loc[, 1:2])
colnames(mesh_points) <- c("x", "y")
ggplot(mesh_points, aes(x = x, y = y)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Easting (km)", y = "Northing (km)", title = "Mesh Locations")
heights <- data$altitude # Altitude of the stations
precipitation <- data$precipitation / 1000 # Precipitation data
data_height_prec <- data.frame(altitude = heights, y = precipitation) # Data frame of altitude and precipitation

station_locations <- cbind(data$x, data$y)
A <- inla.spde.make.A(mesh = mesh, loc = station_locations) # Maps observations to mesh
A2 <- fm_basis(mesh, loc = station_locations) # This does the same
# Checks if dgcmatrices A A2 are equal
all.equal(A, A2)

# Plots precipitation against altitude gg
ggplot(data_height_prec, aes(x = altitude, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)
# Creates a linear model for precipitation against altitude
lm_prec <- lm(y ~ altitude, data = data_height_prec)
beta0 <- lm_prec$coefficients[[1]]
beta1 <- lm_prec$coefficients[[2]]
y <- precipitation - beta0 - beta1 * heights

# Plots the residuals of the linear model
ggplot(data_height_prec, aes(x = altitude, y = y)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE)




# Defines the upper bounds for the quantiles
rho0 <- 10 # Controls the size of kappa
a0 <- 10 # Controls the size of v
sigma_u0 <- 3 # controls standard deviation of field
sigma_epsilon0 <- 3 # control standard deviation of noise
alpha <- 0.05 # Defines the quantile
m_u <- 0 # Setting mean of the field

# Calculates the log prior density function of theta for PC and non-PC priors
log_pc_prior <- log_pc_prior_quantile(
    sigma_u0 = sigma_u0, sigma_epsilon0 = sigma_epsilon0,
    a0 = a0, rho0 = rho0, alpha = alpha
)
# number of iterations
maxit_MAP <- 600

guess <- sim_theta_pc_quantile(
    alpha = alpha, sigma_u0 = sigma_u0,
    sigma_epsilon0 = sigma_epsilon0, a0 = a0, rho0 = rho0, m = 1
)

map_pc <- MAP_prior(
    log_prior = log_pc_prior, mesh = mesh,
    y = y, A = A, m_u = m_u, max_iterations = maxit_MAP,
    theta0 = unlist(guess)
)
map_par <- map_pc$par
v <- c(map_par[2], map_par[3])
v_norm <- sqrt(v[1]^2 + v[2]^2)
v_tilde <- vector_to_half_angle_vector(v)
non_log <- parameters_to_non_log(map_par)
rho <- non_log[[1]]
plt_covariance_of_u(
    log_kappa = map_par[1],
    v = v,
    log_sigma_u = map_par[4],
    l = rho / exp(v_norm / 2), n = 300,
    X_label = "Eastings (km)",
    Y_label = "Northings (km)",
    path = "Manuscript/Simulation_images/covariance_precip.pdf"
)
