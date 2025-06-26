# In this file we plot the anisotropic PC and EG = isotropic PC priors on rho and r = |v| for values from the precipitation study
library(devtools)
library(ggplot2)
library(latex2exp)
library(devtools)
library(tidyr)
library(dplyr)
library(patchwork)
#devtools::document(pkg = "package")
# Plot of prior on kappa
rho0 <- 10 # Controls the size of kappa
a0 <- 10 # Controls the size of v
l <- 300
n_points <- 20000
alpha <- 0.05 # Defines the quantile
path <- "Manuscript/Simulation_images/Priors/Prior_comparison_precip_rho_r.pdf"

plt_pc_priors_kappa_rho_aniso_iso(
  rho0 = rho0,
  a0 = a0,
  alpha = alpha,
  l = l,
  n_points = n_points,
  path = path,
  txt_size = 20
)
