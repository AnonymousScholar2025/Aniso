# library(SPDEaniso)
library(devtools)
library(ggplot2)
library(Matrix)
library(sp)
library(INLA)
library(inlabru)
library(cowplot)
library(purrr)
library(sf)
library(patchwork)
library(latex2exp)
#devtools::document(pkg = "package")


# Hyperparameters on priors -----------------------------------------------

rho0 <- 1 # Controls correlation range
a0 <- 2 # Controls anisotropy
sigma_u0 <- 10 # controls standard deviation of field
sigma_epsilon0 <- 0.5 # control standard deviation of noise

L <- 10 # Size of domain, related to support of range for uniform and beta pripors
width_uniform <- Inf

shape <- 1.1
width_beta <- 20
alpha <- 0.01


# Priors ------------------------------------------------------------------
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
log_uniform_prior <- log_prior_uniform(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  L = L,
  width_support_factor = width_uniform
)


log_beta_prior <- SPDEaniso::log_prior_beta(
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  L = L,
  shape = shape,
  width_support_factor = width_beta
)

# Simulates the "true parameters" from the pc_prior.
set.seed(555)
true_params <- sim_theta_pc_quantile(
  alpha = alpha,
  sigma_u0 = sigma_u0,
  sigma_epsilon0 = sigma_epsilon0,
  a0 = a0,
  rho0 = rho0,
  m = 1
)

# Domain and mesh ---------------------------------------------------------
boundary_sf <- sf::st_sfc(sf::st_polygon(list(rbind(
  c(0, 0.01), c(L, 0.01), c(L, L), c(0, L), c(0, 0.01)
))))
boundary <- fmesher::fm_as_segm(boundary_sf)
mesh <- fmesher::fm_mesh_2d_inla(boundary = boundary, max.edge = c(1, 3))
nodes <- mesh$loc
n <- mesh$n
plot(mesh)
n_observations <- 15
observations <- 10 * matrix(runif(n_observations * 2), ncol = 2)
A <- fmesher::fm_basis(mesh, loc = observations)
m_u <- 0

log_kappa <- true_params$log_kappa
kappa <- exp(log_kappa)
v <- true_params$v
log_sigma_u <- true_params$log_sigma_u
log_sigma_epsilon <- true_params$log_sigma_epsilon
# Sample from noisy data
aniso <- list(rep(kappa, n), matrix(v, n, 2))
x <- fm_aniso_basis_weights_sample(
  x = mesh,
  aniso = aniso,
  log_sigma = log_sigma_u
)
y <- A %*% x + exp(log_sigma_epsilon) * stats::rnorm(nrow(A))

log_priors <- list(
  pc = log_pc_prior,
  EG = log_EG_prior,
  uniform = log_uniform_prior,
  beta = log_beta_prior
)
prior_types <- setNames(as.list(names(log_priors)), names(log_priors))
log_posteriors <- lapply(log_priors, function(log_prior) {
  function(log_kappa,
           v,
           log_sigma_u,
           log_sigma_epsilon) {
    log_posterior_prior(
      log_prior = log_prior,
      mesh = mesh,
      log_kappa = log_kappa,
      v = v,
      log_sigma_epsilon = log_sigma_epsilon,
      log_sigma_u = log_sigma_u,
      y = y,
      A = A,
      m_u = m_u
    )
  }
})


true_par <- c(
  true_params$log_kappa,
  true_params$v,
  true_params$log_sigma_u,
  true_params$log_sigma_epsilon
)
path <- function(i,
                 n_points_log_kappa,
                 n_points_v,
                 l_log_kappa,
                 l_v) {
  paste0(
    "Manuscript/Section_5/Prior_Posterior_Data/pppd_",
    names(log_priors)[i],
    "_n_kp=",
    n_points_log_kappa,
    "_n_v=",
    n_points_v,
    "_l_kp=",
    l_log_kappa,
    "_l_v=",
    l_v,
    ".rds"
  )
}



n_points_log_kappa_pc <- 200
n_points_v_pc <- 100
l_log_kappa_pc <- 5
l_v_pc <- 1
l_log_kappa_uniform <- 5
l_v_uniform <- 5

#Generating data

# df_pc <- prior_posterior_plotter(
#   theta_fixed = true_par,
#   log_priors = log_priors[1],
#   log_posteriors = log_posteriors[1],
#   l_log_kappa = l_log_kappa_pc,
#   l_v = l_v_pc,
#   n_points_log_kappa = n_points_log_kappa_pc,
#   n_points_v = n_points_v_pc,
#   path_data = path(
#     1,
#     n_points_log_kappa_pc,
#     n_points_v_pc,
#     l_log_kappa_pc,
#     l_v_pc
#   ),
#   # path1 = "Manuscript/Simulation_images/prior_posterior_plots_log_kappa.pdf",
#   # path2 = "Manuscript/Simulation_images/prior_posterior_plots_v.pdf"
# )



# df_eg <- prior_posterior_plotter(
#   theta_fixed = true_par,
#   log_priors = log_priors[2],
#   log_posteriors = log_posteriors[2],
#   l_log_kappa = l_log_kappa_pc,
#   l_v = l_v_pc,
#   n_points_log_kappa = n_points_log_kappa_pc,
#   n_points_v = n_points_v_pc,
#   path_data = path(
#     2,
#     n_points_log_kappa_pc,
#     n_points_v_pc,
#     l_log_kappa_pc,
#     l_v_pc
#   )
# )

# df_uniform <- prior_posterior_plotter(
#   theta_fixed = true_par,
#   log_priors = log_priors[3],
#   log_posteriors = log_posteriors[3],
#   l_log_kappa = l_log_kappa_uniform,
#   l_v = l_v_uniform,
#   n_points_log_kappa = n_points_log_kappa_pc,
#   n_points_v = n_points_v_pc,
#   path_data = path(
#     3,
#     n_points_log_kappa_pc,
#     n_points_v_pc,
#     l_log_kappa_uniform,
#     l_v_uniform
#   )
# )


# df_beta <- prior_posterior_plotter(
#   theta_fixed = true_par,
#   log_priors = log_priors[4],
#   log_posteriors = log_posteriors[4],
#   l_log_kappa = l_log_kappa_uniform,
#   l_v = l_v_uniform,
#   n_points_log_kappa = n_points_log_kappa_pc,
#   n_points_v = n_points_v_pc,
#   path_data = path(
#     4,
#     n_points_log_kappa_pc,
#     n_points_v_pc,
#     l_log_kappa_uniform,
#     l_v_uniform
#   )
# )


#Reading data

df_pc <- readRDS(path(
  1,
  n_points_log_kappa_pc,
  n_points_v_pc,
  l_log_kappa_pc,
  l_v_pc
))
df_eg <- readRDS(path(
  2,
  n_points_log_kappa_pc,
  n_points_v_pc,
  l_log_kappa_pc,
  l_v_pc
))
df_uniform <- readRDS(path(
  3,
  n_points_log_kappa_pc,
  n_points_v_pc,
  l_log_kappa_uniform,
  l_v_uniform
))
df_beta <- readRDS(path(
  4,
  n_points_log_kappa_pc,
  n_points_v_pc,
  l_log_kappa_uniform,
  l_v_uniform
))
theta_fixed <- df_pc$theta_fixed
df_list <- list(df_pc, df_eg, df_uniform, df_beta)

path_log_kappa <- paste0(
  "Manuscript/Simulation_images/Priors/prior_posterior_plots_log_kappa_n_log_kappa=",
  n_points_log_kappa_pc,
  "_n_v=",
  n_points_v_pc,
  ".pdf"
)
path_v <- paste0(
  "Manuscript/Simulation_images/Priors/prior_posterior_plots_v_n_log_kappa=",
  n_points_log_kappa_pc,
  "_n_v=",
  n_points_v_pc,
  ".pdf"
)





create_plot <- function(df_v,
                        theta_fixed,
                        prior_types,
                        txt_size = 20,
                        show_colorbar = FALSE) {
  p <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    # geom_vline(
    #   data = data.frame(Parameter = "v1", xintercept = theta_fixed[2]),
    #   aes(xintercept = xintercept),
    #   color = "blue"
    # ) +
    # geom_hline(
    #   data = data.frame(Parameter = "v2", yintercept = theta_fixed[3]),
    #   aes(yintercept = yintercept),
    #   color = "blue"
    # ) +
    facet_wrap(~ FunctionType + PriorType, ncol = length(prior_types)) +
    labs(
      x = expression(v[1]),
      y = expression(v[2]),
      fill = "Density"
    ) +
    theme(text = element_text(size = txt_size), aspect.ratio = 1)

  if (show_colorbar) {
    p <- p + guides(fill = guide_colorbar(barwidth = 20, barheight = 1))
  } else {
    p <- p + guides(fill = "none")
  }

  return(p)
}

plt_kp_v_data_priors_posteriors <- function(df_list,
                                            path_log_kappa = NULL,
                                            path_v = NULL,
                                            txt_size = 12,
                                            width = 10,
                                            height = 5) {
  df_list_1 <- list(df_list[[1]], df_list[[2]])
  df_list_2 <- list(df_list[[3]], df_list[[4]])

  df_log_kappa <- do.call(rbind, lapply(df_list, function(df) {
    df[[1]]
  }))

  df_v_1 <- do.call(rbind, lapply(df_list_1, function(df) {
    df[[2]]
  }))
  df_v_2 <- do.call(rbind, lapply(df_list_2, function(df) {
    df[[2]]
  }))
  p_log_kappa <- ggplot(
    df_log_kappa,
    aes(
      x = x,
      y = value,
      color = PriorType,
      linetype = FunctionType
    )
  ) +
    geom_line(linewidth = 3) +
    # geom_vline(
    #   data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed[1]),
    #   aes(xintercept = xintercept),
    #   color = "blue"
    # ) +
    labs(x = TeX("$\\log(\\kappa)$"), y = "Density") +
    theme(text = element_text(size = txt_size)) +
    theme(legend.position = "bottom") +
    guides(
      linetype = guide_legend(override.aes = list(linewidth = 2, size = 20)),
      color = guide_legend(override.aes = list(size = 4))
    )

  p_v_1 <- create_plot(df_v_1, theta_fixed, prior_types, txt_size, show_colorbar = TRUE)
  p_v_2 <- create_plot(df_v_2, theta_fixed, prior_types, txt_size, show_colorbar = FALSE)

  p_v_combined <- p_v_1 / p_v_2 + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(p_log_kappa)
  print(p_v_combined)
  if (!is.null(path_log_kappa)) {
    ggsave(path_log_kappa, p_log_kappa, width = width, height = height)
  }
  if (!is.null(path_v)) {
    ggsave(path_v, p_v_combined, width = width, height = height)
  }
}

plt_kp_v_data_priors_posteriors(
  df_list,
  path_log_kappa = path_log_kappa,
  path_v = path_v,
  txt_size = 30,
  width = 20,
  height = 15
)

# df_list <- list(df_pc)

# df_log_kappa <- lapply(df_list, function(df) {
#   df[[1]]
# })
# df_log_kappa <- do.call(rbind, df_log_kappa)
# df_v <- lapply(df_list, function(df) {
#   df[[2]]
# })
# df_v <- do.call(rbind, df_v)
# p_log_kappa <- ggplot(df_log_kappa,
#                       aes(
#                         x = x,
#                         y = value,
#                         color = PriorType,
#                         linetype = FunctionType
#                       )) +
#   geom_line() +
#   geom_vline(
#     data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed["log_kappa"]),
#     aes(xintercept = xintercept),
#     color = "blue"
#   ) +
#   labs(x = "log_kappa", y = "Density")

# # Plot for v
# p_v <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
#   geom_tile() +
#   geom_vline(
#     data = data.frame(Parameter = "v1", xintercept = theta_fixed["v1"]),
#     aes(xintercept = xintercept),
#     color = "blue"
#   ) +
#   geom_hline(
#     data = data.frame(Parameter = "v2", yintercept = theta_fixed["v2"]),
#     aes(yintercept = yintercept),
#     color = "blue"
#   ) +
#   facet_wrap( ~ FunctionType + PriorType, ncol = length(prior_types)) +
#   labs(x = "v1", y = "v2", fill = "Density")
# # Print the plots
# print(p_log_kappa)
# print(p_v)
# df_list2 <- list(df_uniform, df_beta)
# df_log_kappa2 <- lapply(df_list2, function(df) {
#   df[[1]]
# })
# df_log_kappa2 <- do.call(rbind, df_log_kappa2)
# df_v2 <- lapply(df_list2, function(df) {
#   df[[2]]
# })
# df_v2 <- do.call(rbind, df_v2)
# p_log_kappa2 <- ggplot(df_log_kappa2,
#                        aes(
#                          x = x,
#                          y = value,
#                          color = PriorType,
#                          linetype = FunctionType
#                        )) +
#   geom_line() +
#   geom_vline(
#     data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed["log_kappa"]),
#     aes(xintercept = xintercept),
#     color = "blue"
#   ) +
#   labs(x = "log_kappa", y = "Density")

# # Plot for v
# p_v2 <- ggplot(df_v2, aes(x = x, y = y, fill = value)) +
#   geom_tile() +
#   geom_vline(
#     data = data.frame(Parameter = "v1", xintercept = theta_fixed["v1"]),
#     aes(xintercept = xintercept),
#     color = "blue"
#   ) +
#   geom_hline(
#     data = data.frame(Parameter = "v2", yintercept = theta_fixed["v2"]),
#     aes(yintercept = yintercept),
#     color = "blue"
#   ) +
#   facet_wrap( ~ FunctionType + PriorType, ncol = length(prior_types)) +
#   labs(x = "v1", y = "v2", fill = "Density")
# # Print the plots
# print(p_log_kappa2)
# print(p_v2)

# p_v / p_v2 + plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")

# p_log_kappa2 + p_log_kappa









# df_list <- list(df_pc, df_uniform, df_beta)
# df_log_kappa <- lapply(df_list, function(df) {
#   df[[1]]
# })
# df_log_kappa <- do.call(rbind, df_log_kappa)
# p_log_kappa <- ggplot(df_log_kappa,
#                        aes(
#                          x = x,
#                          y = value,
#                          color = PriorType,
#                          linetype = FunctionType
#                        )) +
#   geom_line() +
#   geom_vline(
#     data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed["log_kappa"]),
#     aes(xintercept = xintercept),
#     color = "blue"
#   ) +
#   labs(x = "log_kappa", y = "Density")

# p_log_kappa


# #
# #
# # df <- readRDS("Simulation_results/prior_posterior_plots_data_60.rds")
# # true_par <- list(
# #   log_kappa = true_params$log_kappa,
# #   v1 = true_params$v[1],
# #   v2 = true_params$v[2],
# #   log_sigma_u = true_params$log_sigma_u,
# #   log_sigma_epsilon = true_params$log_sigma_epsilon
# # )
# # plot_df_priors_posteriors(df,
# #                           theta_true = true_par,
# #                           path1 = "Manuscript/Simulation_images/prior_posterior_plots_log_kappa.pdf",
# #                           path2 = "Manuscript/Simulation_images/prior_posterior_plots_v.pdf")
# #
# #
# # #
# # # # The plots are different when rho0 and a0 are large. eg (2,4) more for (4,10). But are still quite similar.
# # # n_points_log_kappa <- 10
# # # prior_posterior_plotter_just_log_kappa(
# # #   theta_fixed = map_pc$par,
# # #   log_priors = log_priors,
# # #   log_posteriors = log_posteriors,
# # #   l_log_kappa = -1,
# # #   l_v = 6,
# # #   n_points_log_kappa = n_points_log_kappa,
# # #   n_points_v = 4,
# # #   path_data = paste0(
# # #     "Simulation_results/pppd_log_kappa",
# # #     n_points_log_kappa,
# # #     ".rds"
# # #   ),
# # #   path = paste0(
# # #     "Manuscript/Simulation_images/prior_posterior_plots_just_log_kappa",
# # #     n_points_log_kappa,
# # #     ".pdf"
# # #   )
# # # )
# #
# #
