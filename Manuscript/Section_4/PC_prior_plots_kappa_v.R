# In this file we plot the prior on kappa and v for given values
# of the flexibility parameters lambda, lambda1
library(SPDEaniso)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
if (!requireNamespace("latex2exp", quietly = TRUE)) {
  install.packages("latex2exp")
}
library(latex2exp)
#devtools::document(pkg = "package")
lambda <- 0.1
lambda1 <- 0.1

plt_pc_prior_kappa_v(lambda, lambda1, n_points = 60, path = "Manuscript/Simulation_images/Priors/PC_prior_kappa_v.pdf", txt_size = 20)


PC_prior_log_kappa(1, lambda, lambda1)

curve(PC_prior_log_kappa(x, lambda, lambda1), from = -5, to = 5, n = 1000, xlab = expression(paste("log(", kappa, ")")), ylab = "Density", main = "PC prior on log(kappa)", lwd = 2)
