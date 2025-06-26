# In this file we plot the prior on kappa and v for given values
# of the flexibility parameters lambda, lambda1
library(devtools)
library(ggplot2)
library(latex2exp)
library(devtools)
library(tidyr)
library(dplyr)
library(patchwork)
document()
lambda <- 0.1
lambda1 <- 0.1

plt_pc_prior_kappa_v(lambda, lambda1, n_points = 600, path = "Manuscript/Simulation_images/Priors/PC_prior_kappa_v.pdf", txt_size = 20)


PC_prior_log_kappa(1, lambda, lambda1)

curve(PC_prior_log_kappa(x, lambda, lambda1), from = -5, to = 5, n = 1000, xlab = expression(paste("log(", kappa, ")")), ylab = "Density", main = "PC prior on log(kappa)", lwd = 2)
