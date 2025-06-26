# In this file we plot the covariance function K(x,p),K(x,x) for a nonstationary field under different values of v
# We also plot simulations of the field
# We do this by obtaining the precision associated to the FEM method of the SPDE.
# The covariance can then be used by interpolating the covariance function at the points of interest
# The simulations are obtained using a Cholesky decomposition of the precision matrix
library(SPDEaniso)
library(ggplot2)
library(sf)

# Mesh definition ---------------------------------------------------------
kp <- 1
l <- sqrt(8) / kp * exp(0.5) * 1.25
boundary_sf <- st_sfc(st_polygon(list(rbind(
  c(-l, -l), c(l + 1e-3, -l - 1e-3), c(l, l), c(-l, l), c(-l, -l)
))))
boundary <- fmesher::fm_as_segm(boundary_sf)

mesh <- fmesher::fm_mesh_2d_inla(boundary = boundary, max.edge = c(0.4, 0.8))
plot(mesh)


# Kappa functions ---------------------------------------------------
kappa <- function(x) {
  return(kp)
}
kappa_list <- replicate(4, kappa)

# V functions -------------------------------------------------------------
v_list <- list(function(x) {
  c(-1, 0)
}, function(x) {
  c(0, -1)
}, function(x) {
  c(1, 0)
}, function(x) {
  c(0, 1)
})


# Sigma_u values ---------------------------------------------------
sigma_u_list <- c(1, 1, 1, 1)
txt_size <- 40
# Plotting the field ---------------------------------------------------
titles <- c(bquote(bold(v) == .(c(1, 0))), bquote(bold(v) == .(c(0, 1))), bquote(bold(v) == .(c(-1, 0))), bquote(bold(v) == .(c(0, -1))))
# Should return the same as stationary function
plt_covariance_of_u(
  mesh = mesh,
  point = cbind(0, 0),
  boundary = boundary_sf,
  kappa_list = kappa_list,
  v_list = v_list,
  sigma_u_list = sigma_u_list,
  path = "Manuscript//Simulation_images/Field_plots/covariance_v_rotating_2.pdf",
  titles = sapply(titles, deparse),
  n_col = 2,
  txt_size = txt_size
)

plt_u_sample(
  mesh = mesh,
  boundary = boundary_sf,
  kappa_list = kappa_list,
  v_list = v_list,
  sigma_u_list = sigma_u_list,
  path = "Manuscript//Simulation_images/Field_plots/u_sample_v_rotating.pdf",
  titles = sapply(titles, deparse),
  n_col = 2,
  txt_size = txt_size
)

v_list_non_stationary <- list(function(x) {
  alpha <- 2 * atan2(x[2], x[1])
  norm <- sqrt(sum(x^2))
  norm * c(-sin(alpha), cos(alpha))
}, function(x) {
  c(x[1], x[2])
})


set.seed(123)
u <- plt_u_sample(
  mesh = mesh,
  boundary = boundary_sf,
  kappa_list = kappa_list[1:2],
  v_list = v_list_non_stationary,
  sigma_u_list = sigma_u_list[1:2],
  path = "Manuscript//Simulation_images/Field_plots/u_sample_v_non_stationary.pdf",
  titles = "",
  n_col = 2,
  txt_size = txt_size
)

cov <- plt_covariance_of_u(
  mesh = mesh,
  point = cbind(2, 2),
  boundary = boundary_sf,
  kappa_list = kappa_list[1:2],
  v_list = v_list_non_stationary,
  sigma_u_list = sigma_u_list,
  path = "Manuscript//Simulation_images/Field_plots/covariance_v_non_stationary.pdf",
  n_col = 2,
  txt_size = txt_size
)

plt_marginal_variance_of_u(
  mesh = mesh,
  boundary = boundary_sf,
  kappa_list = kappa_list[1:2],
  v_list = v_list_non_stationary,
  sigma_u_list = sigma_u_list,
  path = "Manuscript//Simulation_images/Field_plots/marginal_variance_v_non_stationary.pdf",
  titles = "",
  n_col = 2,
  txt_size = txt_size
)
