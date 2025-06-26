# In this file we plot the results from the simulation study of different priors for the parameters
# (log(kappa),v,log(sigma_u),log(sigma_epsilon)) of the anisotropic SPDE
# That is we plot the distances to the MAP estimate, the lengths of the credible intervals, the probabilities of the true parameter being in the credible interval, the KL divergence, the complexity of the approximation, the K-S test and the k diagnostics.
library(ggplot2)
library(dplyr)
library(tidyr)
library(devtools)
library(patchwork)
library(latex2exp)
library(xtable)

document()

prior_types <- list(
  pc = "pc",
  EG = "EG",
  uniform = "uniform",
  beta = "beta"
)
approximation_types <- list("Gaussian_median", "importance", "importance_smoothed")
approximation_types_names <- c("GSN", "IMP", "IMPS")
approximation_types <- setNames(approximation_types, approximation_types_names)
approximation_indexes <- c(3)
approximation_types <- approximation_types[approximation_indexes]
n_observations <- 15
number_of_loops <- 600
number_of_weights <- 1000

titles <- list(
  pc = expression(paste(bold(theta)^{
    true
  }, " ~ ", pi[PC])),
  EG = expression(paste(bold(theta)^{
    true
  }, " ~ ", pi[EG])),
  uniform = expression(paste(bold(theta)^{
    true
  }, " ~ ", pi[U])),
  beta = expression(paste(bold(theta)^{
    true
  }, " ~ ", pi[beta]))
)
results_list <- list()
for (prior_type in prior_types) {
  true_parameter_distribution <- prior_type
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

  path_suffix <- paste0(
    "_",
    true_parameter_distribution,
    "_n_obs=",
    n_observations,
    "_n_loops=",
    number_of_loops,
    "_n_weights=",
    number_of_weights,
    "_wu=",
    width_uniform,
    "_wb=",
    width_beta
  )
  results_list[[prior_type]] <- readRDS(paste0(
    "Manuscript/Section_5/Simulation_results/results",
    path_suffix,
    ".rds"
  ))
}
width_uniform <- Inf
width_beta <- 20
path_suffix <- paste0(
  "_n_obs=",
  n_observations,
  "_n_loops=",
  number_of_loops,
  "_n_weights=",
  number_of_weights,
  "_wu=",
  width_uniform,
  "_wb=",
  width_beta
)
variable_indexes <- c(1, 2, 3) # log(kappa), v1,v2 are plotted
# variable_indexes <- c(4, 5) # log_sigma_u, log_sigma_epsilon are plotted
#variable_indexes <- c(1, 2, 3, 4, 5) # all parameters are plotted
labeller <- function(variable, value) {
  latex_labels <- list(
    log_kappa = TeX("$\\log(\\kappa)$"),
    v1 = TeX("$v_1$"),
    v2 = TeX("$v_2$"),
    log_sigma_u = TeX("$\\log(\\sigma_u)$"),
    log_sigma_epsilon = TeX("$\\log(\\sigma_\\epsilon)$")
  )
  value <- lapply(value, function(v) {
    latex_labels[[v]]
  })
  return(value)
}
if (length(variable_indexes) == 2) {
  path_suffix <- paste0("only_sigma", path_suffix)
  labeller <- function(variable, value) {
    latex_labels <- list(
      log_sigma_u = TeX("$\\log(\\sigma_u)$"),
      log_sigma_epsilon = TeX("$\\log(\\sigma_\\epsilon)$")
    )
    value <- lapply(value, function(v) {
      latex_labels[[v]]
    })
    return(value)
  }
}
if (length(variable_indexes) == 5) {
  path_suffix <- paste0("all_par", path_suffix)
}
parameter_names <- rownames(results_list$pc[[1]]$pc$credible_intervals$Gaussian_median)[variable_indexes]
txt_size <- 30
width <- 20
height1 <- 15
height2 <- 7.5
height3 <- 40
txt_angle <- 90

plt_distances_to_MAP(
  results_list = results_list,
  prior_types = prior_types,
  path = paste0("Manuscript/Simulation_images/MAP_loop/Distances_to_map", path_suffix, ".pdf"),
  width = width,
  height = height1,
  titles = titles,
  variable_indexes = variable_indexes,
  txt_size = txt_size,
  txt_angle = txt_angle,
  lower_bound = 10^-2,
  labeller = labeller
)
plt_CI_lengths_and_get_mean_lengths(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = approximation_types,
  parameter_names = parameter_names,
  path = paste0("Manuscript/Simulation_images/MAP_loop/CI_lengths_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height1,
  titles = titles,
  txt_size = txt_size,
  lower_bound = 5 * 10^-1,
  labeller = labeller
)
plt_frequency_true_parameter_in_CI(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = approximation_types,
  parameter_names = parameter_names,
  path = paste0("Manuscript/Simulation_images/MAP_loop/within_CI_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height1,
  titles = titles,
  variable_indexes = variable_indexes,
  txt_size = txt_size
)
KL_approx_types <- list(importance = "importance", smoothed_importance = "smoothed_importance")
KL_approx_types <- setNames(KL_approx_types, c("IMP", "IMPS"))
KL_indexes <- c(2)
KL_approx_types <- KL_approx_types[KL_indexes]

plt_KL_and_get_mean_KL(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = KL_approx_types,
  path = paste0("Manuscript/Simulation_images/MAP_loop/KL_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height2,
  titles = titles,
  txt_size = txt_size
)

plt_probabilities(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = approximation_types,
  parameter_names = parameter_names,
  path = paste0("Manuscript/Simulation_images/MAP_loop/probabilities_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height1,
  titles = titles,
  txt_size = txt_size
)
plt_KS(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = approximation_types,
  parameter_names = parameter_names,
  path1 = paste0("Manuscript/Simulation_images/MAP_loop/KS_distance_smoothed", path_suffix, ".pdf"),
  path2 = paste0("Manuscript/Simulation_images/MAP_loop/KS_pvalue_smoothed", path_suffix, ".pdf"),
  path3 = paste0("Manuscript/Simulation_images/MAP_loop/KS_bridge_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height1,
  titles = titles,
  txt_size = txt_size
)

KS_table(ks_results = read.csv("Manuscript/Section_5/Prior_Posterior_Data/KS_results_1.csv"))
KS_table(ks_results = read.csv("Manuscript/Section_5/Prior_Posterior_Data/KS_results_2.csv"))
KS_table(ks_results = read.csv("Manuscript/Section_5/Prior_Posterior_Data/KS_results_3.csv"))
KS_table(ks_results = read.csv("Manuscript/Section_5/Prior_Posterior_Data/KS_results_4.csv"))




plt_complexity_and_get_mean_complexity(
  results_list = results_list,
  prior_types = prior_types,
  approximation_types = approximation_types,
  path = paste0("Manuscript/Simulation_images/MAP_loop/complexity_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height2,
  titles = titles,
  txt_size = txt_size
)
plt_complexity_and_get_mean_complexity(
  results_list = results_list,
  prior_types = prior_types[1:2],
  approximation_types = approximation_types,
  path = paste0("Manuscript/Simulation_images/MAP_loop/complexity_pc_EG_smoothed", path_suffix, ".pdf"),
  width = width,
  height = height2,
  titles = titles,
  txt_size = txt_size
)

plt_k_diagnostics(
  results_list = results_list,
  prior_types = prior_types,
  path = paste0("Manuscript/Simulation_images/MAP_loop/k_diagnostics", path_suffix, ".pdf"),
  width = width,
  height = height2,
  titles = titles
)
