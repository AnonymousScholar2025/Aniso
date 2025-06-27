# In this file we plot the covariance function for a stationary field under different values of v
# Using a transformation of the Matern covariance function
library(SPDEaniso)
library(pracma)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(grid)
log_kappa_list <- rep(0, 4)
v_list <- list(c(1, 0), c(0, 1), c(-1, 0), c(0, -1))
log_sigma_u_list <- rep(0, 4)
n_col <- 4
path <- paste0("Manuscript//Simulation_images/Field_plots/covariance_v_rotating.pdf")
titles <- c(bquote(bold(v) == .(c(1, 0))), bquote(bold(v) == .(c(0, 1))), bquote(bold(v) == .(c(-1, 0))), bquote(bold(v) == .(c(0, -1))))
titles <- sapply(titles, deparse)
plt_covariance_of_u_stationary(
  log_kappa_list = log_kappa_list,
  v_list = v_list,
  log_sigma_u_list = log_sigma_u_list,
  l = sqrt(8) / 1.25 / exp(log_kappa_list[1]) * exp(sqrt(sum(v_list[[1]]^2))),
  n = 300,
  path = path,
  titles,
  n_col = 2,
  txt_size = 20
)
