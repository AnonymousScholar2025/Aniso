# DESCRIPTION: This file contains functions to plot the different priors and posteriors and tests used in the simulation study comparing different priors (distances to MAP, the credible interval lengths. the mean lengths of the credible intervals, the KS test, model complexity, etc.)
# Additionally, there are functions for plotting anisotropic fields and their covariances over a spatial domain.

#' @title Plot the pc prior for r and rho
#' @description Plot the pc prior for r and rho
#' @param rho0 An unexpectedly small correlation range
#' @param a0 An unexpectedly small anisotropy
#' @param alpha The quantile
#' @param l The length of the partition
#' @param n_points The number of points in the partition
#' @param path A character string of the path to save the plot
#' @param txtsize An integer of the text size by default 20
plt_pc_priors_kappa_rho_aniso_iso <- function(rho0,
                                              a0,
                                              alpha = 0.01,
                                              l = 20,
                                              n_points = 1000,
                                              path = NULL,
                                              txt_size = 20) {
  pc_rho_iso <- pc_prior_rho_iso(rho0 = rho0, alpha = alpha)
  pc_prior_rho_aniso <- pc_prior_rho(rho0 = rho0, a0 = a0, alpha = alpha)
  PC_prior_r <- PC_prior_r_quantile(a0 = a0, alpha = alpha)
  EG_prior_r <- EG_prior_r(a0 = a0, alpha = alpha)

  priors_rho <- list(anisotropic = pc_prior_rho_aniso, isotropic = pc_prior_rho_iso)
  priors_r <- list(PC_prior_r = PC_prior_r, EG_prior_r = EG_prior_r)
  partition_rho <- seq(0.01, l, length.out = n_points)
  partition_r <- seq(0.01, log(a0) * 2, length.out = n_points)


  df_rho <- data.frame(
    partition_rho = rep(partition_rho, 2),
    value = c(
      pc_rho_iso(partition_rho),
      pc_prior_rho_aniso(partition_rho)
    ),
    prior_type = rep(c("EG", "PC"), each = length(partition_rho))
  )

  df_r <- data.frame(
    partition_r = rep(partition_r, 2),
    value = c(PC_prior_r(partition_r), EG_prior_r(partition_r)),
    prior_type = rep(c("PC", "EG"), each = length(partition_r))
  )

  df_rho_long <- df_rho %>%
    pivot_longer(
      cols = -c(partition_rho, prior_type),
      names_to = "variable",
      values_to = "value"
    )
  df_r_long <- df_r %>%
    pivot_longer(
      cols = -c(partition_r, prior_type),
      names_to = "variable",
      values_to = "value"
    )



  p_rho <- ggplot(
    df_rho_long,
    aes(x = partition_rho, y = value, color = prior_type)
  ) +
    geom_line() +
    scale_color_manual(values = c("blue", "red")) +
    labs(
      x = TeX("$\\rho$"),
      y = "Density",
      color = "Prior Type"
    )
  p_r <- ggplot(df_r_long, aes(x = partition_r, y = value, color = prior_type)) +
    geom_line() +
    scale_color_manual(values = c("blue", "red")) +
    labs(
      x = TeX("$r$"),
      y = "Density",
      color = "Prior Type"
    )

  p <- combined_plot <- wrap_plots(list(p_rho, p_r))
  p <- p + plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      text = element_text(size = txt_size)
    )
  print(p)
  if (!is.null(path)) {
    ggsave(
      path,
      plot = p,
      width = 10,
      height = 5,
      dpi = 300
    )
  }
}
#' @title Plot the pc priors on kappa and v
#' @description Plot the pc priors on kappa and v
#' @param lambda Controls the size of kappa
#' @param lambda1 Controls the size of v
#' @param alpha The quantile
#' @param l The length of the partition
#' @param n_points The number of points in the partition
#' @param path A character string of the path to save the plot
#' @param txt_size An integer of the text size by default 20

plt_pc_prior_kappa_v <- function(lambda,
                                 lambda1,
                                 alpha = 0.01,
                                 l = 4,
                                 n_points = 100,
                                 path = NULL,
                                 txt_size = 20) {
  pc_prior_kappa <- PC_prior_kappa(lambda = lambda, lambda1 = lambda1)
  pc_prior_v <- PC_prior_v(lambda1 = lambda1)

  partition_kappa <- seq(0.01, l, length.out = n_points)
  partition_v <- expand.grid(
    v1 = seq(-l, l, length.out = n_points),
    v2 = seq(-l, l, length.out = n_points)
  )

  df_kappa <- data.frame(
    partition_kappa = partition_kappa,
    value = pc_prior_kappa(partition_kappa)
  )

  df_v <- data.frame(
    partition_v1 = partition_v$v1,
    partition_v2 = partition_v$v2,
    value = apply(partition_v, 1, pc_prior_v)
  )

  df_kappa_long <- df_kappa %>%
    pivot_longer(
      cols = -c(partition_kappa),
      names_to = "variable",
      values_to = "value"
    )
  df_v_long <- df_v %>%
    pivot_longer(
      cols = -c(partition_v1, partition_v2),
      names_to = "variable",
      values_to = "value"
    )

  p_kappa <- ggplot(df_kappa_long, aes(x = partition_kappa, y = value)) +
    geom_line() +
    labs(
      x = TeX(r'($\kappa$)'),
      y = TeX(r'($\pi_{\kappa}(\kappa)$)')
    ) +
    theme(
      text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0), units = "inches")
    )
  p_v <- ggplot(df_v_long, aes(x = partition_v1, y = partition_v2, group = interaction(partition_v1, partition_v2))) +
    geom_tile(aes(fill = value)) +
    geom_contour(aes(z = value), color = "black") +
    scale_fill_gradientn(colors = c("white", "red")) +
    labs(
      x = TeX("$v_1$"),
      y = TeX("$v_2$"),
      fill = "Density"
    ) +
    theme(
      legend.key.height = unit(1.5, "cm"),
      legend.text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0.5), units = "inches")
    ) +
    coord_equal()

  p <- wrap_plots(list(p_kappa, p_v))
  p <- p + plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      text = element_text(size = txt_size)
    )
  print(p)
  if (!is.null(path)) {
    ggsave(
      path,
      plot = p,
      width = 10.25,
      height = 4,
      dpi = 300
    )
  }
}
# Covariance of solution to (kappa^2-nabla\cdot H_v\nabla)u=\sqrt{4\pi}\sigma_u W_noise
# Anisotropic matern with marginal variance sigma_u
plt_covariance_of_u_stationary <- function(log_kappa_list,
                                           v_list,
                                           log_sigma_u_list,
                                           l = 4,
                                           n = 300,
                                           path = NULL,
                                           titles = NULL,
                                           X_label = expression(x[1]),
                                           Y_label = expression(x[2]),
                                           legend_title = "K",
                                           X_breaks = c(-4, 0, 4),
                                           Y_breaks = c(-4, 0, 4),
                                           n_col = 4,
                                           n_row = 2,
                                           legend_height = 4,
                                           txt_size = 30) {
  # Define Matern covariance function
  plots <- list()
  all_data <- data.frame()
  for (i in seq_along(log_kappa_list)) {
    log_kappa <- log_kappa_list[[i]]
    v <- v_list[[i]]
    sigma_u <- exp(log_sigma_u_list[[i]])

    kappa <- exp(log_kappa)

    correlation_matern <- function(x_norm) {
      kappa_x_norm <- kappa * x_norm
      sigma_u^2 * (kappa_x_norm) * besselK(kappa_x_norm, 1)
    }
    v_norm <- sqrt(v[1]^2 + v[2]^2)
    # Normalised half angle vector
    if (v_norm == 0) {
      vhn <- c(1, 0)
    } else {
      vhn <- vector_to_half_angle_vector(v) / v_norm
    }

    vhn_o <- c(-vhn[2], vhn[1])

    # Define covariance of u, K(x,y)=r(|x-y|_{F^{-1}}) where F= sqrt(H)/kappa has orthonormal
    # eigenvectors {vhn,vhn_o} and eigenvalues {exp(-|v|/2),exp(-|v|/2)}
    K_kappa_v <- function(x) {
      deformed_norm_x <- sqrt(exp(-v_norm) * (x %*% vhn)^2 + exp(v_norm) * (x %*% vhn_o)^
        2)
      correlation_matern(deformed_norm_x)
    }

    # Stores data in a data frame
    pxl <- expand.grid(
      x = seq(-l, l, length.out = n),
      y = seq(-l, l, length.out = n)
    )
    pxl$K <- mapply(function(x, y) {
      K_kappa_v(c(x, y))
    }, pxl$x, pxl$y)
    pxl$plot <- i # add a new variable indicating the plot number
    all_data <- rbind(all_data, pxl)
  }
  # labels <- setNames(titles, seq_along(titles))
  # labeller_func <- as_labeller(labels, default = label_parsed)
  p <- ggplot(all_data) +
    geom_tile(aes(x = x, y = y, fill = K)) +
    scale_x_continuous(breaks = X_breaks) +
    scale_y_continuous(breaks = Y_breaks) +
    geom_contour(aes(x = x, y = y, z = K)) +
    scale_fill_gradientn(colors = c("white", "red")) +
    labs(x = X_label, y = Y_label, fill = legend_title) +
    facet_wrap(~plot, ncol = n_col, as.table = FALSE) + # create facets
    coord_equal() +
    theme(
      text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
      strip.text = element_blank(),
      legend.key.height = unit(0.5 * legend_height / n_col, "cm")
    )

  print(p)

  if (!is.null(path)) {
    ggsave(path, height = 4 * n_row * n_col, dpi = 300)
  }
}

#' @ title Plot realizations of a field on a mesh
#' @description Plot realizations of a list of fields on a mesh given a list of functions for kappa and v
#' @param mesh A mesh object
#' @param boundary Boundary of the mesh
#' @param kappa_list A list of functions for kappa
#' @param v_list A list of functions for v
#' @param sigma_u_list A list of values for the standard deviation sigma_u
#' @param path A character string of the path to save the plot
#' @param titles A character vector of titles for each plot

plt_u_sample <- function(mesh,
                         boundary,
                         kappa_list,
                         v_list,
                         sigma_u_list,
                         path = NULL,
                         titles,
                         n_col = 4,
                         txt_size = 30) {
  nodes <- mesh$loc
  all_data <- data.frame()
  # labels <- setNames(titles, seq_along(titles))
  # labeller_func <- as_labeller(labels, default = label_parsed)
  for (i in seq_along(kappa_list)) {
    kappa <- kappa_list[[i]]
    v <- v_list[[i]]
    sigma_u <- sigma_u_list[i]
    kappa_values <- apply(nodes, 1, kappa)
    vec_values <- t(apply(nodes, 1, v))
    aniso <- list(kappa_values, vec_values)
    u_mesh <- sigma_u * fm_aniso_sample(mesh, aniso) # Sample of field on mesh
    pxl <- fm_pixels(mesh, dims = c(150, 150), mask = boundary) # Define pixels for plotting masking outside the boundary
    pxl$u <- fm_evaluate(mesh, # Evaluate the field at the pixels and store in the data frame
      loc = pxl, field = u_mesh
    )
    pxl$plot <- i
    all_data <- rbind(all_data, pxl)
  }
  all_data$x_coord <- sapply(all_data$geometry, function(g) {
    g[1]
  })
  all_data$y_coord <- sapply(all_data$geometry, function(g) {
    g[2]
  })
  p <- ggplot(all_data) +
    geom_tile(aes(geometry = geometry, fill = u), stat = "sf_coordinates") +
    scale_fill_gradientn(
      colours = c("#0000FFFF", "#FFFFFFFF", "#FF0000FF"),
      limits = c(-max(abs(all_data$u)), max(abs(all_data$u)))
    ) +
    scale_x_continuous(breaks = c(-4, 0, 4)) +
    scale_y_continuous(breaks = c(-4, 0, 4)) +
    facet_wrap(~plot, ncol = n_col, as.table = FALSE) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    coord_equal() +
    theme(
      text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
      strip.text = element_blank(),
      legend.key.height = unit(1.5, "cm")
    )
  if (!is.null(path)) {
    ggsave(path,
      height = 5 * 4 / n_col,
      width = 10.5,
      dpi = 300
    )
  }
  p
}
#' @title Plot the Covariance of the field with a point
#' @description Plot the covariance K(x,p) of the field with a point for a list of functions for kappa and v
#' @param mesh A mesh object
#' @param boundary Boundary of the mesh
#' @param point A point to evaluate the covariance with
#' @param kappa_list A list of values for kappa (not its logarithm)
#' @param v_list A list of vectors for v
#' @param sigma_u_list A list of values for the standard deviation sigma_u
#' @param path A character string of the path to save the plot
#' @param titles A character vector of titles for each plot
#' @param n_col An integer of the number of columns for the facet wrap

plt_covariance_of_u <- function(mesh,
                                boundary = NULL,
                                point = cbind(0, 0),
                                kappa_list,
                                v_list,
                                sigma_u_list,
                                path = NULL,
                                titles = NULL,
                                n_col = 4,
                                txt_size = 30) {
  all_data <- data.frame()
  nodes <- mesh$loc
  for (i in seq_along(kappa_list)) {
    kappa <- kappa_list[[i]]
    v <- v_list[[i]]
    sigma_u <- sigma_u_list[i]
    kappa_values <- apply(nodes, 1, kappa)
    vec_values <- t(apply(nodes, 1, v))
    aniso <- list(kappa_values, vec_values)
    Q <- fm_aniso_precision(mesh, aniso)
    A1 <- Matrix::Diagonal(nrow(Q)) # Identity matrix to evaluate the field at the mesh. A1 u = u
    A2 <- fm_basis(mesh, point) # Matrix (vector) to evaluate the field at p A2u = u(p)
    # Calculates the covariance of A1u with A2u where u has precision Q. That is the covariance of u(x) with u(p)
    K <- sigma_u^2 * fm_covariance(Q, A1, A2)
    pxl <- fm_pixels(mesh, dims = c(150, 150), mask = boundary)
    pxl$K <- fm_evaluate(mesh, loc = pxl, field = as.vector(K))
    pxl$plot <- i
    all_data <- rbind(all_data, pxl)
  }
  all_data$X <- sf::st_coordinates(all_data)[, 1]
  all_data$Y <- sf::st_coordinates(all_data)[, 2]
  p <- ggplot(all_data) +
    geom_tile(aes(X, Y, fill = K), alpha = 1) +
    geom_contour(aes(X, Y, z = K)) +
    scale_fill_gradientn(
      colours = c("#FFFFFFFF", "#FF0000FF"),
      limits = c(min(all_data$K), max(abs(all_data$K)))
    ) +
    scale_x_continuous(breaks = c(-4, 0, 4)) +
    scale_y_continuous(breaks = c(-4, 0, 4)) +
    facet_wrap(~plot, ncol = n_col, as.table = FALSE) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    coord_equal() +
    theme(
      text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
      strip.text = element_blank(),
      legend.key.height = unit(1 * 4 / n_col, "cm")
    )

  if (!is.null(path)) {
    ggsave(path,
      height = 5 * 4 / n_col,
      width = 10.25,
      dpi = 300
    )
  }
  p
}

#' @title Plot the marginal variance K(x,x) of the field
#' @description Plot the covariance K(x,x) of the field with a point for a list of functions for kappa and v
#' @param mesh A mesh object
#' @param boundary Boundary of the mesh
#' @param kappa_list A list of values for kappa
#' @param v_list A list of vectors for v
#' @param sigma_u_list A list of values for the standard deviation sigma_u
#' @param path A character string of the path to save the plot
#' @param titles A character vector of titles for each plot
#' @param n_col An integer of the number of columns for the facet wrap

plt_marginal_variance_of_u <- function(mesh,
                                       boundary,
                                       kappa_list,
                                       v_list,
                                       sigma_u_list,
                                       path = NULL,
                                       titles,
                                       n_col = 4,
                                       txt_size = 30) {
  all_data <- data.frame()
  labels <- setNames(titles, seq_along(titles))
  labeller_func <- as_labeller(labels, default = label_parsed)
  nodes <- mesh$loc
  for (i in seq_along(kappa_list)) {
    kappa <- kappa_list[[i]]
    v <- v_list[[i]]
    sigma_u <- sigma_u_list[i]
    kappa_values <- apply(nodes, 1, kappa)
    vec_values <- t(apply(nodes, 1, v))
    aniso <- list(kappa_values, vec_values)
    Q <- fm_aniso_precision(mesh, aniso)
    V <- diag(INLA::inla.qinv(Q), nrow(Q), ncol(Q))
    pxl <- fm_pixels(mesh, dims = c(150, 150), mask = boundary)
    pxl$V <- fm_evaluate(mesh, loc = pxl, field = as.vector(sqrt(V)))^
      2
    pxl$plot <- i
    all_data <- rbind(all_data, pxl)
  }
  p <- ggplot(all_data) +
    geom_tile(aes(geometry = geometry, fill = V),
      stat = "sf_coordinates",
      alpha = 1
    ) +
    scale_fill_gradientn(
      colours = c("#FFFFFFFF", "#FF0000FF"),
      limits = c(0, max(abs(all_data$V)))
    ) +
    facet_wrap(~plot,
      ncol = n_col,
      labeller = labeller_func,
      as.table = FALSE
    ) +
    scale_x_continuous(breaks = c(-4, 0, 4)) +
    scale_y_continuous(breaks = c(-4, 0, 4)) +
    facet_wrap(~plot, ncol = n_col, as.table = FALSE) +
    labs(x = expression(x[1]), y = expression(x[2])) +
    coord_equal() +
    theme(
      text = element_text(size = txt_size),
      plot.margin = unit(c(0, 0, 0, 0), units = "inches"),
      strip.text = element_blank(),
      legend.key.height = unit(1 * 4 / n_col, "cm")
    )
  print(p)
  if (!is.null(path)) {
    ggsave(path,
      height = 5 * 4 / n_col,
      width = 10.25,
      dpi = 300
    )
  }
}




# Define a labelling function
latex_labeller <- function(variable, value) {
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

#' @title Plot distances to MAP
#' @description Plot distances to MAP for each prior type
#' @param results_list A list of results from the simulation containing the distances to MAP
#' @param prior_types A character vector of prior types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot. By default, it is 600
#' @param width An integer of the width of the plot. By default, it is 10
#' @param height An integer of the height of the plot. By default, it is 10
#' @param titles A character vector of titles for each plot
#' @param variable_indexes An integer vector of the indexes of the variables to plot. By default, it is c(1, 2, 3, 4, 5)
#' @param txt_size An integer of the text size. By default, it is 20
#' @return A ggplot object
#' @export

plt_distances_to_MAP <- function(results_list,
                                 prior_types,
                                 path = NULL,
                                 dpi = 100,
                                 width = 10,
                                 height = 5,
                                 titles,
                                 labeller = latex_labeller,
                                 variable_indexes = c(1, 2, 3, 4, 5),
                                 txt_size = 20,
                                 lower_bound = 10^-2,
                                 txt_angle = 45) {
  # Define variable names in order
  variable_names <- c("log_kappa", "v1", "v2", "log_sigma_u", "log_sigma_epsilon")

  # Subset variable names based on indexes
  selected_variables <- variable_names[variable_indexes]

  # Initialize a list to store the plots
  plots <- list()

  # Loop over the list of results
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    all_distances <- data.frame()
    for (prior_type in prior_types) {
      distances_to_MAP <- lapply(results, function(x) {
        x[[prior_type]]$distance_vector[variable_indexes]
      })
      distances_to_MAP <- do.call(rbind, distances_to_MAP)
      distances_to_MAP <- as.data.frame(distances_to_MAP)
      distances_to_MAP$iteration <- seq_len(nrow(distances_to_MAP))
      distances_to_MAP <- reshape2::melt(distances_to_MAP, id.vars = "iteration")
      distances_to_MAP$prior_type <- prior_type
      all_distances <- rbind(all_distances, distances_to_MAP)
    }

    # Update the 'variable' column with the selected variable names
    all_distances$variable <- factor(all_distances$variable, labels = selected_variables)

    p <- ggplot(all_distances) +
      stat_ecdf(aes(value, color = prior_type, linetype = prior_type)) +
      labs(color = "Prior", linetype = "Prior") +
      scale_x_continuous(trans = "log10") +
      coord_cartesian(xlim = c(lower_bound, NA)) +
      facet_wrap(~variable, labeller = labeller) +
      labs(title = titles[[i]]) +
      labs(x = "Distance to MAP", y = "CDF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size),
        axis.text.x = element_text(angle = txt_angle)
      )
    if (i %in% c(2, 4)) {
      p <- p + ylab(NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    if (i %in% c(1, 2)) {
      p <- p + xlab(NULL)
    }

    plots[[i]] <- p
  }

  # Combine the plots
  combined_plot <- wrap_plots(plots)

  # Show the legend only once and place it at the bottom
  combined_plot <- combined_plot + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  height <- if (length(variable_indexes) > 3) height else 0.7 * height

  if (!is.null(path)) {
    ggsave(
      path,
      plot = combined_plot,
      dpi = dpi,
      height = height,
      width = width
    )
  }

  return(combined_plot)
}


#' @title Mean distances to MAP and standard deviation of Gaussian approximation around the median
#' @description Gets the mean distances to MAP and the standard deviation of the Gaussian approximation around the median for each prior type
#' @param results A list of results from the simulation containing the distances to MAP
#' @param prior_types A character vector of prior types
#' @return A list of vectors containing the mean distances to MAP and the standard deviation of the Gaussian approximation around the median for each prior type
#' @export

mean_distance_to_MAP_and_std_dev_of_Gaussian_approximation <- function(results, prior_types) {
  mean_distances <- lapply(prior_types, function(prior_type) {
    mean_distances <- sapply(1:5, function(i) {
      all_distances <- sapply(seq_along(results), function(j) {
        results[[j]][[prior_type]]$distance_vector[i]
      })
      mean(all_distances)
    })
    names(mean_distances) <- parameter_names

    std_dev_estimates_Gaussian_median <- do.call(
      rbind,
      lapply(results, function(x) {
        x[[prior_type]]$std_dev_estimates_Gaussian_median
      })
    )
    mean_std_dev <- colMeans(std_dev_estimates_Gaussian_median)
    list(mean_distances = mean_distances, mean_std_dev = mean_std_dev)
  })
  mean_distances
}


#' @title Plot credible interval lengths and get mean lengths
#' @description Plot credible interval lengths and get mean lengths for each prior type and approximation type
#' @param results A list of results from the simulation containing the credible intervals
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @param titles A character vector of titles for each plot
#' @param txt_size An integer of the text size by default 20
#' @return A list of data frames containing the mean lengths for each prior type and approximation type
#' @export



plt_CI_lengths_and_get_mean_lengths <- function(results_list,
                                                prior_types,
                                                approximation_types,
                                                parameter_names,
                                                path = NULL,
                                                dpi = 100,
                                                width = 10,
                                                height = 10,
                                                titles,
                                                labeller = latex_labeller,
                                                txt_size = 20,
                                                lower_bound = 5 * 10^-1) {
  # Initialize a list to store the plots
  plots <- list()

  # Initialize a list to store the mean lengths
  mean_lengths_list <- list()

  # Loop over the list of results
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    lengths_df <- lapply(prior_types, function(prior_type) {
      lapply(approximation_types, function(approximation_type) {
        lengths <- lapply(parameter_names, function(parameter_name) {
          all_lengths <- sapply(seq_along(results), function(j) {
            length <- diff(results[[j]][[prior_type]]$credible_intervals[[approximation_type]][parameter_name, ])
            length[length == 0] <- 10^-3 # Remove 0 lengths for the log transformation
            length
          })
          all_lengths
        })
        lengths <- do.call(cbind, lengths)
        colnames(lengths) <- parameter_names
        lengths
      })
    })
    mean_lengths_df <- lapply(lengths_df, function(prior_type) {
      lapply(prior_type, function(approximation_type) {
        colMeans(approximation_type)
      })
    })
    all_lengths <- data.frame()
    for (prior_type in prior_types) {
      for (approximation_type in names(approximation_types)) {
        lengths <- lengths_df[[prior_type]][[approximation_type]]
        lengths <- as.data.frame(lengths)
        lengths$iteration <- seq_len(nrow(lengths))
        lengths <- reshape2::melt(lengths, id.vars = "iteration") # Necessary to use ggplot as it expects a data frame in long format
        lengths$prior_type <- prior_type
        lengths$approximation_type <- approximation_type
        all_lengths <- rbind(all_lengths, lengths)
      }
    }

    p <- ggplot(all_lengths) +
      stat_ecdf(aes(value, color = prior_type, linetype = prior_type)) +
      labs(color = "Prior", linetype = "Prior") +
      facet_wrap(~variable, labeller = labeller) +
      scale_x_continuous(trans = "log10", limits = c(0.5, NA)) +
      coord_cartesian(xlim = c(lower_bound, NA)) +
      labs(x = "Credible Interval Length", y = "CDF") +
      labs(title = titles[[i]]) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )
    if (i %in% c(2, 4)) {
      p <- p + ylab(NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    if (i %in% c(1, 2)) {
      p <- p + xlab(NULL)
    }

    # Add the plot to the list
    plots[[i]] <- p

    # Add the mean lengths to the list
    mean_lengths_list[[prior_types[[i]]]] <- mean_lengths_df
  }

  # Combine the plots
  combined_plot <- wrap_plots(plots, ncol = length(plots))

  # Show the legend only once and place it at the bottom
  combined_plot <- combined_plot + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(combined_plot)

  height <- if (length(parameter_names) > 3) height else 0.7 * height

  # Save the combined plot
  if (!is.null(path)) {
    ggsave(
      path,
      plot = combined_plot,
      dpi = dpi,
      width = width,
      height = height
    )
  }

  # Return the combined plot and the list of mean lengths
  mean_lengths_list
}


#' @title Plot the frequency of the true parameter being in the credible interval
#' @description Plot the frequency of the true parameter being in the credible interval for each prior type and approximation type
#' @param results A list of results from the simulation containing the credible intervals
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @param titles A character vector of titles for each plot
#' @param txt_size An integer of the text size by default 20
#' @export

# Define a mapping function
approximation_type_mapper <- function(approximation_type) {
  approximation_type_abbreviations <- list(
    Gaussian_median = "GSN",
    importance = "IMP",
    importance_smoothed = "IMPS"
  )
  return(approximation_type_abbreviations[[approximation_type]])
}

plt_frequency_true_parameter_in_CI <- function(results_list,
                                               prior_types = prior_types,
                                               approximation_types = approximation_types,
                                               parameter_names = parameter_names,
                                               path = NULL,
                                               dpi = 100,
                                               width = 10,
                                               height = 10,
                                               titles,
                                               variable_indexes = c(1, 2, 3, 4, 5),
                                               txt_size = 20) {
  # Initialize a list to store the plots
  plots <- list()
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    within_ci <- lapply(prior_types, function(prior_type) {
      lapply(approximation_types, function(approximation_type) {
        within_binary <- sapply(results, function(x) {
          x[[prior_type]][["true_parameter_within_c_interval"]][[approximation_type]][variable_indexes]
        })
        list(
          mean = rowMeans(within_binary),
          variance = apply(within_binary, 1, var) / length(results)
        )
      })
    })

    # Convert the list to a data frame
    within_ci <- do.call(rbind, lapply(names(within_ci), function(prior_type) {
      do.call(rbind, lapply(names(within_ci[[prior_type]]), function(approximation_type) {
        data.frame(
          prior_type = prior_type,
          approximation_type = approximation_type,
          # Use the mapping function
          parameter = factor(
            names(within_ci[[prior_type]][[approximation_type]]$mean),
            levels = c(
              "log_kappa",
              "v1",
              "v2",
              "log_sigma_u",
              "log_sigma_epsilon"
            )
          ),
          mean = c(within_ci[[prior_type]][[approximation_type]]$mean),
          std_dev = sqrt(c(within_ci[[prior_type]][[approximation_type]]$variance))
        )
      }))
    }))
    # Create a named vector of labels
    labels <- c(
      log_kappa = TeX("$\\log(\\kappa)$"),
      v1 = TeX("$v_1$"),
      v2 = TeX("$v_2$"),
      log_sigma_u = TeX("$\\log(\\sigma_u)$"),
      log_sigma_epsilon = TeX("$\\log(\\sigma_\\epsilon)$")
    )

    # Use the labels in your plot
    p <- ggplot(
      within_ci,
      aes(
        x = parameter,
        y = mean,
        ymin = mean - 2 * std_dev,
        ymax = mean + 2 * std_dev,
        color = prior_type,
        linetype = prior_type
      )
    ) +
      geom_pointrange() +
      geom_text(aes(label = round(mean, 2)), vjust = 0, hjust = -.25, size = 10) +
      labs(x = "Parameter", y = expression(paste("Frequency of ", theta^
        true, " in CI"))) +
      labs(title = titles[[i]]) +
      theme(
        axis.text.x = element_text(
          size = txt_size
        ),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      ) +
      scale_y_continuous(trans = "log10") +
      scale_x_discrete(labels = labels) # Use the labels
    if (i %in% c(2, 4)) {
      p <- p + ylab(NULL) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    if (i %in% c(1, 2)) {
      p <- p + xlab(NULL)
    }

    plots[[i]] <- p
  }
  # Combine the plots
  combined_plot <- wrap_plots(plots)

  # Show the legend only once and place it at the bottom
  combined_plot <- combined_plot + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  # Save the combined plot
  if (!is.null(path)) {
    ggsave(
      path,
      plot = combined_plot,
      dpi = dpi,
      width = width,
      height = height
    )
  }

  # Return the combined plot
  return(combined_plot)
}

#' @title Plot KL divergences
#' @description Plot KL divergences for each prior type and approximation type
#' @param results A list of results from the simulation containing the KL divergences
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @param titles A character vector of titles for each plot
#' @param txt_size An integer of the text size by default 20
#' @export

plt_KL_and_get_mean_KL <- function(results_list,
                                   prior_types,
                                   approximation_types,
                                   path = NULL,
                                   dpi = 100,
                                   width = 10,
                                   height = 10,
                                   titles,
                                   txt_size = 20) {
  library(ggplot2)
  library(patchwork)
  library(reshape2)

  plots <- list()
  mean_KL_list <- list()
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    KL <- lapply(prior_types, function(prior_type) {
      lapply(approximation_types, function(approximation_type) {
        sapply(results, function(x) {
          kl_values <- x[[prior_type]]$importance[[paste0(
            "KL_divergence_",
            approximation_type,
            "_Gaussian_median"
          )]]
          kl_values
        })
      })
    })
    KL <- lapply(KL, function(prior_type) {
      lapply(prior_type, function(approximation_type) {
        max_value <- if (any(is.finite(unlist(approximation_type)))) {
          max(unlist(approximation_type)[is.finite(unlist(approximation_type))])
        } else {
          10^4
        }
        replace(
          approximation_type,
          is.na(unlist(approximation_type)) |
            unlist(approximation_type) == Inf,
          max_value
        )
      })
    })

    all_KL <- data.frame()
    for (prior_type in prior_types) {
      for (approximation_type in names(approximation_types)) {
        KL_divergences <- KL[[prior_type]][[approximation_type]]
        KL_divergences <- as.data.frame(KL_divergences)
        KL_divergences$iteration <- seq_len(nrow(KL_divergences))
        KL_divergences <- reshape2::melt(KL_divergences, id.vars = "iteration")
        KL_divergences$prior_type <- prior_type
        KL_divergences$approximation_type <- approximation_type
        all_KL <- rbind(all_KL, KL_divergences)
      }
    }

    p <- ggplot(all_KL) +
      stat_ecdf(aes(value, color = prior_type, linetype = prior_type)) +
      labs(color = "Prior", linetype = "Prior") +
      scale_x_continuous(trans = "log10") +
      labs(title = titles[[i]]) +
      labs(x = "KL Divergence", y = "CDF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )

    if (i != 1) {
      p <- p + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    }

    plots[[i]] <- p

    KL_Gaussian_median_mean <- lapply(prior_types, function(prior_type) {
      lapply(names(approximation_types), function(approximation_type) {
        mean(KL[[prior_type]][[approximation_type]])
      })
    })
    mean_KL_list[[i]] <- KL_Gaussian_median_mean
  }

  combined_plot <- wrap_plots(plots, ncol = 4) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = txt_size),
      axis.title.x = element_text(size = txt_size)
    )

  print(combined_plot)
  height <- if (length(parameter_names) > 3) height else 0.7 * height
  if (!is.null(path)) {
    ggsave(
      path,
      plot = combined_plot,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  KL_Gaussian_median_mean
}


#' @title Plot probabilities
#' @description Plot probabilities for each prior type and approximation type that the marginal posterior is smaller than the true parameter value
#' @param results A list of results from the simulation containing the probabilities
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_probabilities <- function(results_list,
                              prior_types,
                              approximation_types,
                              parameter_names,
                              path = NULL,
                              dpi = 100,
                              width = 10,
                              height = 10,
                              titles,
                              labeller = latex_labeller,
                              txt_size = 20) {
  plots <- list()
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    all_probabilities <- data.frame()
    for (prior_type in prior_types) {
      for (approximation_key in names(approximation_types)) {
        approximation_type <- approximation_types[[approximation_key]]
        probabilities <- sapply(results, function(x) {
          x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]]
        })
        for (j in seq_along(parameter_names)) {
          # Use the names of approximation types as labels
          df <- data.frame(
            prob = unlist(probabilities[j, ]),
            parameter = parameter_names[[j]],
            prior = prior_type,
            approximation = approximation_key
          )
          all_probabilities <- rbind(all_probabilities, df)
        }
      }
    }

    all_probabilities$parameter <- factor(all_probabilities$parameter, levels = parameter_names)
    p <- ggplot(all_probabilities) +
      stat_ecdf(aes(prob, color = prior, linetype = prior)) +
      labs(color = "Prior", linetype = "Prior") +
      geom_abline(
        slope = 1,
        intercept = 0,
        color = "black"
      ) +
      facet_wrap(~parameter, labeller = labeller) +
      labs(title = titles[[i]]) +
      labs(x = expression(u[i]), y = "CDF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )
    if (i %in% c(2, 4)) {
      p <- p + ylab(NULL)
    }

    if (i %in% c(1, 2)) {
      p <- p + xlab(NULL)
    }
    if (i %in% c(2, 4)) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    plots[[i]] <- p
  }
  combined_plot <- wrap_plots(plots)
  combined_plot <- combined_plot + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(combined_plot)
  height <- if (length(parameter_names) > 3) height else 0.7 * height
  if (!is.null(path)) {
    ggsave(path,
      dpi = dpi,
      width = width,
      height = height
    )
  }
}

#' @title Plot KS statistics and p-values
#' @description Plot KS statistics and p-values for each prior type and approximation type
#' @param results A list of results from the simulation containing the probabilities
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param parameter_names A character vector of parameter names
#' @param path1 A character string of the path to save the plot for the KS statistics
#' @param path2 A character string of the path to save the plot for the p-values
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

inla_KS_copy <- function(x, y, ...) {
  if (any(is.na(x))) {
    x <- x[!is.na(x)]
  }
  test <- ks.test(x, y, ...)
  n <- length(x)
  Fn <- ((1:n) - 0.5) / n
  CDF_true <- y(sort(x))
  empirical.diff <- (Fn - CDF_true) * sqrt(n)
  # T <- max(abs(empirical.diff))
  # ylim <- c(-1, 1) * max(1, T)
  plot(
    CDF_true,
    empirical.diff,
    type = "l",
    xlim = c(0, 1),
    ylim = ylim,
    main = paste("K-S-test, p-value = ", test$p.value),
    ylab = "(Fn-F) sqrt(n)",
    xlab = "Quantile"
  )
  # lines
  lines(c(0, 1), c(0, 0), type = "l")
  # Showing ellipse
  lines(Fn, 2 * sqrt(Fn * (1 - Fn)), type = "l")
  lines(Fn, -2 * sqrt(Fn * (1 - Fn)), type = "l")
  # lines(c(0, 1, NA, 0, 1), c(T, T, NA, -T, -T), type = "l")
  invisible(test)
}

latex_labeller2 <- function(variable, value) {
  latex_labels <- c(
    log_kappa = expression(paste(log(kappa))),
    v1 = expression(paste(v[1])),
    v2 = expression(paste(v[2])),
    log_sigma_u = expression(paste(log(sigma[u]))),
    log_sigma_epsilon = expression(paste(log(sigma[epsilon])))
  )
  return(latex_labels)
}
plt_KS <- function(results_list,
                   prior_types,
                   approximation_types,
                   parameter_names,
                   path1 = NULL,
                   path2 = NULL,
                   path3 = NULL,
                   dpi = 100,
                   width = 10,
                   height = 10,
                   titles,
                   txt_size = 20) {
  plots1 <- list()
  plots2 <- list()
  plots3 <- list()
  for (j in seq_along(results_list)) {
    results <- results_list[[j]]
    all_probabilities <- data.frame()
    for (prior_type in prior_types) {
      for (approximation_key in names(approximation_types)) {
        approximation_type <- approximation_types[[approximation_key]]
        probabilities <- sapply(results, function(x) {
          x[[prior_type]][["importance"]][[paste0("probabilities_", approximation_type)]]
        })
        for (i in seq_along(parameter_names)) {
          df <- data.frame(
            prob = unlist(probabilities[i, ]),
            parameter = parameter_names[[i]],
            prior = prior_type,
            approximation = approximation_key
          )
          all_probabilities <- rbind(all_probabilities, df)
        }
      }
    }
    KS_results <- data.frame()
    n <- length(results)
    Fn <- ((1:n) - 0.5) / n
    lines_data <- data.frame(
      Fn = Fn,
      upper = 2 * sqrt(Fn * (1 - Fn)),
      lower = -2 * sqrt(Fn * (1 - Fn))
    )
    for (i in seq_along(parameter_names)) {
      # Get the probabilities for the current parameter
      probabilities <- all_probabilities[all_probabilities$parameter == parameter_names[[i]], ]
      # Calculate the KS statistic for each prior and approximation type
      for (prior_type in prior_types) {
        for (approximation_key in names(approximation_types)) {
          x <- probabilities[probabilities$prior == prior_type &
            probabilities$approximation == approximation_key, ]$prob
          y <- punif
          CDF_true <- y(sort(x))
          empirical_diff <- (Fn - CDF_true) * sqrt(n)
          KS_result <- ks.test(x, y)
          # Add the KS statistic and p-value for the current parameter to the data frame
          KS_results <- rbind(KS_results, suppressWarnings(
            data.frame(
              parameter = parameter_names[[i]],
              prior = prior_type,
              approximation = approximation_key,
              statistic = KS_result$statistic,
              p_value = KS_result$p.value,
              KS_x = CDF_true,
              KS_y = empirical_diff
            )
          ))
        }
      }
    }
    KS_results$parameter <- factor(KS_results$parameter, levels = parameter_names)
    labels <- c(
      log_kappa = TeX("$\\log(\\kappa)$"),
      v1 = TeX("$v_1$"),
      v2 = TeX("$v_2$"),
      log_sigma_u = TeX("$\\log(\\sigma_u)$"),
      log_sigma_epsilon = TeX("$\\log(\\sigma_\\epsilon)$")
    )
    write.csv(KS_results, file = paste0("Manuscript/Section_5/Prior_Posterior_Data/KS_results_", j, ".csv"), row.names = FALSE)

    p1 <- ggplot(KS_results) +
      geom_point(aes(
        x = parameter,
        y = statistic,
        color = prior,
        shape = prior
      )) +
      labs(color = "Prior", shape = "Prior") +
      scale_y_continuous(trans = "log10") +
      scale_x_discrete(labels = labels) +
      labs(title = titles[[j]]) +
      labs(x = "Parameter", y = "KS Statistic") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )
    if (j %in% c(2, 4)) {
      p1 <- p1 + ylab(NULL)
    }
    if (j %in% c(1, 2)) {
      p1 <- p1 + xlab(NULL)
    }
    plots1[[j]] <- p1


    p2 <- ggplot(KS_results) +
      geom_point(aes(
        x = parameter,
        y = p_value,
        color = prior,
        shape = prior
      )) +
      labs(color = "Prior", shape = "Prior") +
      scale_y_continuous(trans = "log10") +
      scale_x_discrete(labels = labels) +
      labs(title = titles[[j]]) +
      labs(x = "Parameter", y = "p-value") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )
    if (j %in% c(2, 4)) {
      p2 <- p2 + ylab(NULL)
    }
    if (j %in% c(1, 2)) {
      p2 <- p2 + xlab(NULL)
    }

    plots2[[j]] <- p2
    KS_results$parameter <- factor(all_probabilities$parameter, levels = parameter_names)
    p3 <- ggplot(KS_results) +
      geom_line(aes(
        x = KS_x,
        y = KS_y,
        color = prior,
        linetype = prior
      )) +
      labs(color = "Prior", linetype = "Prior") +
      geom_line(
        data = lines_data,
        aes(x = Fn, y = upper),
        color = "black"
      ) +
      geom_line(
        data = lines_data,
        aes(x = Fn, y = lower),
        color = "black"
      ) +
      facet_wrap(~parameter, labeller = latex_labeller) +
      labs(title = titles[[j]]) +
      labs(x = "Quantile", y = "(Fn-F) sqrt(n)") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )
    if (j %in% c(2, 4)) {
      p3 <- p3 + ylab(NULL)
    }

    if (j %in% c(1, 2)) {
      p3 <- p3 + xlab(NULL)
      p3 <- p3 + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    plots3[[j]] <- p3
  }
  combined_plot1 <- wrap_plots(plots1)
  combined_plot2 <- wrap_plots(plots2)
  combined_plot3 <- wrap_plots(plots3)
  combined_plot1 <- combined_plot1 + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  combined_plot2 <- combined_plot2 + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  combined_plot3 <- combined_plot3 + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  print(combined_plot1)
  height <- if (length(parameter_names) > 3) height else 0.7 * height

  if (!is.null(path1)) {
    ggsave(path1,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  print(combined_plot2)
  if (!is.null(path2)) {
    ggsave(path2,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  print(combined_plot3)
  if (!is.null(path3)) {
    ggsave(path3,
      dpi = dpi,
      width = width,
      height = height
    )
  }

  # # We only get the ones using last approximation type
  # KS_results <- KS_results[KS_results$approximation == approximation_types[length(approximation_types)], ]
  # # Since we only have one approximation type, we can eliminate it from the names
  # KS_results$approximation <- NULL
  # # aa <- KS_results %>%
  # #     mutate(result = paste0("statistic: ", round(statistic, 2), ", p_value: ", round(p_value, 2)))

  # # # Convert the table from long format to wide format
  # # wide_table <- aa %>%
  # #     pivot_wider(names_from = "parameter", values_from = "result")


  # # # Convert the wide table to a LaTeX table
  # # latex_table <- xtable(KS_results)

  # # # Print the LaTeX table
  # # print(latex_table, type = "latex", include.rownames = FALSE)
  # # KS_results
}

#' @title Plot complexity and get mean complexity
#' @description Plot complexity and get mean complexity for each prior type and approximation type
#' @param results A list of results from the simulation containing the complexity
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @return A list of vectors containing the mean complexity for each prior type and approximation type
#' @export

plt_complexity_and_get_mean_complexity <- function(results_list,
                                                   prior_types,
                                                   approximation_types,
                                                   path = NULL,
                                                   dpi = 100,
                                                   width = 40,
                                                   height = 10,
                                                   titles,
                                                   txt_size = 20) {
  plots <- list()
  mean_complexity_list <- list()
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    complexity <- lapply(prior_types, function(prior_type) {
      lapply(approximation_types, function(approximation_type) {
        sapply(results, function(x) {
          complexity <- x[[prior_type]]$importance[[paste0("complexity_", approximation_type)]]
          complexity
        })
      })
    })

    all_complexity <- data.frame()
    for (prior_type in prior_types) {
      for (approximation_type in names(approximation_types)) {
        complexity_values <- complexity[[prior_type]][[approximation_type]]
        complexity_values <- as.data.frame(complexity_values)
        complexity_values$iteration <- seq_len(nrow(complexity_values))
        complexity_values <- reshape2::melt(complexity_values, id.vars = "iteration")
        complexity_values$prior_type <- prior_type
        complexity_values$approximation_type <- approximation_type
        all_complexity <- rbind(all_complexity, complexity_values)
      }
    }

    p <- ggplot(all_complexity) +
      stat_ecdf(aes(value, color = prior_type, linetype = prior_type)) +
      labs(color = "Prior", linetype = "Prior") +
      scale_x_continuous(trans = "log10") +
      labs(title = titles[[i]]) +
      labs(x = "Complexity", y = "CDF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = txt_size)
      )

    # Remove y-axis labels and ticks for all plots except the first one
    if (i != 1) {
      p <- p + theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    }

    plots[[i]] <- p

    complexity_mean <- lapply(prior_types, function(prior_type) {
      lapply(names(approximation_types[2:3]), function(approximation_type) {
        mean(complexity[[prior_type]][[approximation_type]])
      })
    })
  }

  combined_plot <- wrap_plots(plots, ncol = length(results_list)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  print(combined_plot)
  height <- if (length(parameter_names) > 3) height else 0.7 * height
  if (!is.null(path)) {
    ggsave(
      path,
      plot = combined_plot,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  complexity_mean
}



#' @title Plot k diagnostics
#' @description Plot k diagnostics for each prior type
#' @param results A list of results from the simulation containing the k diagnostics
#' @param prior_types A character vector of prior types
#' @param approximation_types A character vector of approximation types
#' @param path A character string of the path to save the plot
#' @param dpi An integer of the dpi to save the plot
#' @param width An integer of the width of the plot
#' @param height An integer of the height of the plot
#' @export

plt_k_diagnostics <- function(results_list,
                              prior_types,
                              path = NULL,
                              dpi = 100,
                              width = 10,
                              height = 10,
                              titles) {
  plots <- list()
  for (i in seq_along(results_list)) {
    results <- results_list[[i]]
    all_k_diagnostics <- data.frame(matrix(ncol = length(prior_types), nrow = length(results)))
    colnames(all_k_diagnostics) <- prior_types
    for (prior_type in prior_types) {
      all_k_diagnostics[[prior_type]] <- sapply(results, function(x) {
        x[[prior_type]]$importance$k_diagnostic
      })
    }
    all_k_diagnostics <- as.data.frame(all_k_diagnostics)
    all_k_diagnostics$iteration <- seq_len(nrow(all_k_diagnostics))
    all_k_diagnostics <- reshape2::melt(all_k_diagnostics, id.vars = "iteration")

    p <- ggplot(all_k_diagnostics) +
      stat_ecdf(aes(value, color = variable)) +
      labs(color = "Prior", color = "Variable") +
      labs(title = titles[[i]]) +
      labs(x = "k diagnostic", y = "CDF") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)
      )

    if (i %in% c(2, 4)) {
      p <- p + ylab(NULL)
    }

    if (i %in% c(1, 2)) {
      p <- p + xlab(NULL)
    }
    if (i %in% c(2, 4)) {
      p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }

    plots[[i]] <- p
  }
  combined_plot <- wrap_plots(plots)
  combined_plot <- combined_plot + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  if (!is.null(path)) {
    ggsave(path,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  print(combined_plot)
}

plt_weights_cdf <- function(results,
                            prior_types,
                            path = NULL,
                            dpi = 100,
                            width = 10,
                            height = 10) {
  all_weights <- list()
  for (prior_type in prior_types) {
    weights_unsmoothed <- c(unlist(
      sapply(results, function(x) {
        x[[prior_type]]$importance$log_unnormalized_weights
      })
    ))
    weights_smoothed <- c(unlist(
      sapply(results, function(x) {
        x[[prior_type]]$importance$log_unnormalized_weights_smoothed
      })
    ))
    weights <- data.frame(
      prior_type = prior_type,
      weights_unsmoothed = weights_unsmoothed,
      weights_smoothed = weights_smoothed
    )
    all_weights[[prior_type]] <- weights
  }
  all_weights_df <- do.call(rbind, all_weights)

  # Reshape the data to long format
  all_weights_long <- pivot_longer(
    all_weights_df,
    c(weights_unsmoothed, weights_smoothed),
    names_to = "weight_type",
    values_to = "weight"
  )

  p <- ggplot(all_weights_long) +
    stat_ecdf(aes(weight, color = prior_type, linetype = weight_type)) +
    labs(x = "Log weight", y = "Cumulative density") +
    theme(legend.position = "bottom") +
    xlim(c(-5, 0))

  if (!is.null(path)) {
    ggsave(path,
      dpi = dpi,
      width = width,
      height = height
    )
  }
  print(p)
}

# # Function for plotting marginal prior and posterior densities for log_kappa
# prior_posterior_plotter_just_log_kappa <- function(theta_fixed = map_pc$par, log_priors = log_priors,
#                                                    log_posteriors = log_posteriors, l_log_kappa = 3, l_v = 6, n_points_log_kappa = 10, n_points_v = 5, path_data = NULL,
#                                                    path = NULL) {
#   make_evaluable_prior_and_posterior <- function(log_prior) {
#     function(log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon) {
#       tryCatch(
#         {
#           log_prior(
#             log_kappa = log_kappa, v = c(v1, v2),
#             log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
#           )
#         },
#         error = function(e) {
#           # Return -Inf when an error occurs
#           -Inf
#         }
#       )
#     }
#   }

#   log_priors_2 <- lapply(log_priors, make_evaluable_prior_and_posterior)
#   log_posteriors_2 <- lapply(log_posteriors, make_evaluable_prior_and_posterior)
#   priors_and_posteriors <- list(prior = log_priors_2, posterior = log_posteriors_2)


#   partitions <- lapply(1:3, function(i) {
#     if (i == 1) {
#       seq(theta_fixed[i] - l_log_kappa, theta_fixed[i] + l_log_kappa, length.out = n_points_log_kappa)
#     } else {
#       seq(theta_fixed[i] - l_v, theta_fixed[i] + l_v, length.out = n_points_v)
#     }
#   })

#   grid <- expand.grid(
#     log_kappa = partitions[[1]], v1 = partitions[[2]], v2 = partitions[[3]],
#     log_sigma_u = theta_fixed[4], log_sigma_epsilon = theta_fixed[5]
#   )

#   df_list_log_kappa <- list()
#   for (function_type in names(priors_and_posteriors)) {
#     for (prior_type in names(priors_and_posteriors[[function_type]])) {
#       print(paste("Evaluating", function_type, prior_type))
#       pi <- priors_and_posteriors[[function_type]][[prior_type]]
#       grid$value <- pmap_dbl(grid, ~
#         pi(..1, ..2, ..3, ..4, ..5),
#       .progress = TRUE
#       )
#       marginal_log_kappa <- aggregate(value ~ log_kappa, data = grid, FUN = function(x) sum(exp(x)))
#       # Normalize to 1 as max
#       marginal_log_kappa$value <- marginal_log_kappa$value / max(marginal_log_kappa$value)

#       # Store the results in the lists
#       df_list_log_kappa[[paste(function_type, prior_type)]] <- data.frame(
#         x = marginal_log_kappa$log_kappa,
#         value = marginal_log_kappa$value,
#         FunctionType = function_type,
#         PriorType = prior_type,
#         stringsAsFactors = FALSE
#       )
#     }
#   }

#   # Combine the lists of data frames into data frame
#   df_log_kappa <- do.call(rbind, df_list_log_kappa)

#   # Plot for log_kappa
#   p_log_kappa <- ggplot(df_log_kappa, aes(x = x, y = value, color = PriorType, linetype = FunctionType)) +
#     geom_line() +
#     geom_vline(
#       data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed["log_kappa"]),
#       aes(xintercept = xintercept), color = "blue"
#     ) +
#     labs(x = "log_kappa", y = "Density")

#   print(p_log_kappa)

#   if (!is.null(path_data)) {
#     saveRDS(list(df_log_kappa, theta_fixed), path_data)
#   }
#   if (!is.null(path)) {
#     ggsave(path, p_log_kappa, device = "pdf", dpi = 100, width = 10, height = 5)
#   }
# }

# prior_posterior_plotter_just_v <- function(theta_fixed = map_pc$par, log_priors = log_priors,
#                                            log_posteriors = log_posteriors, l_log_kappa = 3, l_v = 6, n_points_log_kappa = 10, n_points_v = 10, path_data = NULL,
#                                            path = NULL) {
#   make_evaluable_prior_and_posterior <- function(log_prior) {
#     function(log_kappa, v1, v2, log_sigma_u, log_sigma_epsilon) {
#       tryCatch(
#         {
#           log_prior(
#             log_kappa = log_kappa, v = c(v1, v2),
#             log_sigma_u = log_sigma_u, log_sigma_epsilon = log_sigma_epsilon
#           )
#         },
#         error = function(e) {
#           # Return -Inf when an error occurs
#           -Inf
#         }
#       )
#     }
#   }

#   log_priors_2 <- lapply(log_priors, make_evaluable_prior_and_posterior)
#   log_posteriors_2 <- lapply(log_posteriors, make_evaluable_prior_and_posterior)
#   priors_and_posteriors <- list(prior = log_priors_2, posterior = log_posteriors_2)


#   partitions <- lapply(1:3, function(i) {
#     if (i == 1) {
#       seq(theta_fixed[i] - l_log_kappa, theta_fixed[i] + l_log_kappa, length.out = n_points)
#     } else {
#       seq(theta_fixed[i] - l_v, theta_fixed[i] + l_v, length.out = n_points)
#     }
#   })

#   grid <- expand.grid(
#     log_kappa = partitions[[1]], v1 = partitions[[2]], v2 = partitions[[3]],
#     log_sigma_u = theta_fixed[4], log_sigma_epsilon = theta_fixed[5]
#   )

#   df_list_v <- list()
#   start_time <- Sys.time()

#   for (function_type in names(priors_and_posteriors)) {
#     for (prior_type in names(priors_and_posteriors[[function_type]])) {
#       print(paste("Evaluating", function_type, prior_type))
#       pi <- priors_and_posteriors[[function_type]][[prior_type]]
#       grid$value <- pmap_dbl(grid, ~
#         pi(..1, ..2, ..3, ..4, ..5),
#       .progress = TRUE
#       )
#       marginal_v <- aggregate(value ~ v1 + v2, data = grid, FUN = function(x) sum(exp(x)))
#       # Normalize to 1 as max
#       marginal_v$value <- marginal_v$value / max(marginal_v$value)

#       # Store the results in the lists
#       df_list_v[[paste(function_type, prior_type)]] <- data.frame(
#         x = marginal_v$v1,
#         y = marginal_v$v2,
#         value = marginal_v$value,
#         FunctionType = function_type,
#         PriorType = prior_type,
#         stringsAsFactors = FALSE
#       )
#     }
#   }

#   # Combine the lists of data frames into two data frames
#   df_v <- do.call(rbind, df_list_v)

#   # Plot for log_kappa

#   # Plot for v
#   p_v <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
#     geom_tile() +
#     geom_vline(
#       data = data.frame(Parameter = "v1", xintercept = theta_fixed["v1"]),
#       aes(xintercept = xintercept), color = "blue"
#     ) +
#     geom_hline(
#       data = data.frame(Parameter = "v2", yintercept = theta_fixed["v2"]),
#       aes(yintercept = yintercept), color = "blue"
#     ) +
#     facet_wrap(~ FunctionType + PriorType, ncol = length(prior_types)) +
#     labs(x = "v1", y = "v2", fill = "Density")
#   # Print the plots
#   print(p_v)

#   if (!is.null(path_data)) {
#     saveRDS(list(df_log_kappa, df_v, theta_fixed), path_data)
#   }
#   if (!is.null(path)) {
#     ggsave(path2, p_v, device = "pdf", dpi = 100, width = 10, height = 5)
#   }
# }



prior_posterior_plotter <- function(theta_fixed = map_pc$par,
                                    log_priors = log_priors,
                                    log_posteriors = log_posteriors,
                                    l_log_kappa = 3,
                                    l_v = 6,
                                    n_points_log_kappa = 100,
                                    n_points_v = 10,
                                    path_data = NULL,
                                    path1 = NULL,
                                    path2 = NULL) {
  make_evaluable_prior_and_posterior <- function(log_prior) {
    function(log_kappa,
             v1,
             v2,
             log_sigma_u,
             log_sigma_epsilon) {
      tryCatch(
        {
          log_prior(
            log_kappa = log_kappa,
            v = c(v1, v2),
            log_sigma_u = log_sigma_u,
            log_sigma_epsilon = log_sigma_epsilon
          )
        },
        error = function(e) {
          # Return -Inf when an error occurs-Inf
        }
      )
    }
  }

  log_priors_2 <- lapply(log_priors, make_evaluable_prior_and_posterior)
  log_posteriors_2 <- lapply(log_posteriors, make_evaluable_prior_and_posterior)
  priors_and_posteriors <- list(prior = log_priors_2, posterior = log_posteriors_2)


  partitions <- lapply(1:3, function(i) {
    if (i == 1) {
      seq(theta_fixed[i] - l_log_kappa,
        theta_fixed[i] + l_log_kappa,
        length.out = n_points_log_kappa
      )
    } else {
      seq(theta_fixed[i] - l_v, theta_fixed[i] + l_v, length.out = n_points_v)
    }
  })

  grid <- expand.grid(
    log_kappa = partitions[[1]],
    v1 = partitions[[2]],
    v2 = partitions[[3]],
    log_sigma_u = theta_fixed[4],
    log_sigma_epsilon = theta_fixed[5]
  )

  df_list_log_kappa <- list()
  df_list_v <- list()
  counter <- 0
  start_time <- Sys.time()

  for (function_type in names(priors_and_posteriors)) {
    for (prior_type in names(priors_and_posteriors[[function_type]])) {
      print(paste("Evaluating", function_type, prior_type))
      pi <- priors_and_posteriors[[function_type]][[prior_type]]
      grid$value <- pmap_dbl(grid, ~
        pi(..1, ..2, ..3, ..4, ..5), .progress = TRUE)
      marginal_log_kappa <- aggregate(
        value ~ log_kappa,
        data = grid,
        FUN = function(x) {
          sum(exp(x))
        }
      )
      marginal_v <- aggregate(
        value ~ v1 + v2,
        data = grid,
        FUN = function(x) {
          sum(exp(x))
        }
      )
      # replace NA with 0
      marginal_log_kappa$value[is.na(marginal_log_kappa$value)] <- 0
      marginal_v$value[is.na(marginal_v$value)] <- 0
      # Normalize to 1 as max
      marginal_log_kappa$value <- marginal_log_kappa$value / max(marginal_log_kappa$value)
      marginal_v$value <- marginal_v$value / max(marginal_v$value)

      # Store the results in the lists
      df_list_log_kappa[[paste(function_type, prior_type)]] <- data.frame(
        x = marginal_log_kappa$log_kappa,
        value = marginal_log_kappa$value,
        FunctionType = function_type,
        PriorType = prior_type,
        stringsAsFactors = FALSE
      )
      df_list_v[[paste(function_type, prior_type)]] <- data.frame(
        x = marginal_v$v1,
        y = marginal_v$v2,
        value = marginal_v$value,
        FunctionType = function_type,
        PriorType = prior_type,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine the lists of data frames into two data frames
  df_log_kappa <- do.call(rbind, df_list_log_kappa)
  df_v <- do.call(rbind, df_list_v)

  # Plot for log_kappa
  p_log_kappa <- ggplot(
    df_log_kappa,
    aes(
      x = x,
      y = value,
      color = PriorType,
      linetype = FunctionType
    )
  ) +
    geom_line() +
    geom_vline(
      data = data.frame(Parameter = "log_kappa", xintercept = theta_fixed["log_kappa"]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    labs(x = "log_kappa", y = "Density")

  # Plot for v
  p_v <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    geom_vline(
      data = data.frame(Parameter = "v1", xintercept = theta_fixed["v1"]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    geom_hline(
      data = data.frame(Parameter = "v2", yintercept = theta_fixed["v2"]),
      aes(yintercept = yintercept),
      color = "blue"
    ) +
    facet_wrap(~ FunctionType + PriorType, ncol = length(prior_types)) +
    labs(x = "v1", y = "v2", fill = "Density")
  # Print the plots
  print(p_log_kappa)
  print(p_v)
  results <- list(
    log_kappa = df_log_kappa,
    v = df_v,
    theta_fixed = theta_fixed
  )

  if (!is.null(path_data)) {
    saveRDS(results, path_data)
  }
  if (!is.null(path1)) {
    ggsave(
      path1,
      p_v,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }
  if (!is.null(path2)) {
    ggsave(
      path2,
      p_v,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }
  return(results)
}

plot_df_priors_posteriors <- function(df,
                                      path1 = NULL,
                                      path2 = NULL,
                                      theta_true = map_pc$par) {
  df_log_kappa <- df[[1]]
  df_v <- df[[2]]
  p_log_kappa <- ggplot(
    df_log_kappa,
    aes(
      x = x,
      y = value,
      color = PriorType,
      linetype = FunctionType
    )
  ) +
    geom_line() +
    geom_vline(
      data = data.frame(Parameter = "log_kappa", xintercept = theta_true[[1]]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    labs(x = "log_kappa", y = "Density")

  p_v <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    geom_vline(
      data = data.frame(Parameter = "v1", xintercept = theta_true[[2]]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    geom_hline(
      data = data.frame(Parameter = "v2", yintercept = theta_true[[3]]),
      aes(yintercept = yintercept),
      color = "blue"
    ) +
    facet_wrap(~ FunctionType + PriorType, ncol = 4) +
    labs(x = "v1", y = "v2", fill = "Density")

  if (!is.null(path1)) {
    ggsave(
      path1,
      p_log_kappa,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }
  if (!is.null(path2)) {
    ggsave(
      path2,
      p_v,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }

  print(p_log_kappa)
  print(p_v)
}

plot_data_prior_posterior <- function(df_list,
                                      path1 = NULL,
                                      path2 = NULL,
                                      theta_true = map_pc$par) {
  df_log_kappa <- df_list[[1]]
  df_v <- df_list[[2]]
  p_log_kappa <- ggplot(
    df_log_kappa,
    aes(
      x = x,
      y = value,
      color = PriorType,
      linetype = FunctionType
    )
  ) +
    geom_line() +
    geom_vline(
      data = data.frame(Parameter = "log_kappa", xintercept = theta_true[[1]]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    labs(x = "log_kappa", y = "Density")

  p_v <- ggplot(df_v, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    geom_vline(
      data = data.frame(Parameter = "v1", xintercept = theta_true[[2]]),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    geom_hline(
      data = data.frame(Parameter = "v2", yintercept = theta_true[[3]]),
      aes(yintercept = yintercept),
      color = "blue"
    ) +
    facet_wrap(~ FunctionType + PriorType, ncol = 4) +
    labs(x = "v1", y = "v2", fill = "Density")

  if (!is.null(path1)) {
    ggsave(
      path1,
      p_log_kappa,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }
  if (!is.null(path2)) {
    ggsave(
      path2,
      p_v,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }

  print(p_log_kappa)
  print(p_v)
}


prior_posterior_plotter_old <- function(theta_fixed = map_pc$par,
                                        log_priors = log_priors,
                                        log_posteriors = log_posteriors,
                                        l = 2,
                                        n_points = 10,
                                        n_parameters_to_plot = 3,
                                        path = NULL) {
  function_types <- list(
    prior = "prior",
    posterior = "posterior",
    path = NULL
  )
  ########## UNNORMALIZED Gaussian_median APPROXIMATION TO THE POSTERIOR############

  ### UNNORMALIZED LOG FUNCTION SO THEY ALL START AT 0###
  unnormalize_prior_and_posterior <- function(log_prior) {
    function(log_kappa,
             v1,
             v2,
             log_sigma_u,
             log_sigma_epsilon) {
      log_prior(
        log_kappa = log_kappa,
        v = c(v1, v2),
        log_sigma_u = log_sigma_u,
        log_sigma_epsilon = log_sigma_epsilon
      ) - log_prior(
        log_kappa = theta_fixed[1],
        v = theta_fixed[2:3],
        log_sigma_u = theta_fixed[4],
        log_sigma_epsilon = theta_fixed[5]
      )
    }
  }

  unnormalized_priors <- lapply(log_priors, unnormalize_prior_and_posterior)
  unnormalized_posteriors <- lapply(log_posteriors, unnormalize_prior_and_posterior)
  unnormalized_priors_and_posteriors <- list(prior = unnormalized_priors, posterior = unnormalized_posteriors)

  # Restricting the functions to one parameter
  restricting_function_to_one_parameter <- function(f, x0) {
    f_list <- lapply(seq_along(x0)[1:n_parameters_to_plot], function(i) {
      function(x) {
        x0_copy <- x0
        x0_copy[i] <- x
        unname(do.call(f, as.list(x0_copy)))
      }
    })
    names(f_list) <- names(x0[1:n_parameters_to_plot])
    f_list
  }
  restricted_priors_and_posteriors <- lapply(unnormalized_priors_and_posteriors, function(f) {
    lapply(f, restricting_function_to_one_parameter, theta_fixed)
  })


  # Getting data for plotting
  partitions <- lapply(seq_along(theta_fixed)[1:n_parameters_to_plot], function(i) {
    seq(theta_fixed[i] - l, theta_fixed[i] + l, length.out = n_points)
  })
  names(partitions) <- names(theta_fixed[1:n_parameters_to_plot])

  plot_data <- do.call(rbind, lapply(names(restricted_priors_and_posteriors), function(function_type) {
    do.call(rbind, lapply(names(restricted_priors_and_posteriors[[function_type]]), function(prior_type) {
      do.call(rbind, lapply(names(
        restricted_priors_and_posteriors[[function_type]][[prior_type]]
      ), function(parameter_name) {
        # Calculate the function values
        x_values <- partitions[[parameter_name]]
        y_values <- sapply(x_values, restricted_priors_and_posteriors[[function_type]][[prior_type]][[parameter_name]])
        # Normalize y values so max is 0
        y_values <- y_values - max(y_values)

        data.frame(
          x = x_values,
          Value = y_values,
          Parameter = parameter_name,
          FunctionType = function_type,
          PriorType = prior_type,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))


  # Create a single ggplot object instead of a list of plots
  p <- ggplot(
    plot_data,
    aes(
      x = x,
      y = exp(Value),
      color = PriorType,
      linetype = FunctionType
    )
  ) +
    geom_line() +
    geom_vline(
      data = data.frame(
        Parameter = names(theta_fixed[1:n_parameters_to_plot]),
        xintercept = unlist(theta_fixed[1:n_parameters_to_plot])
      ),
      aes(xintercept = xintercept),
      color = "blue"
    ) +
    facet_wrap(~Parameter, labeller = latex_labeller, ncol = 2) +
    labs(y = "Density")

  # If you want to save the plot
  if (!is.null(path)) {
    ggsave(
      path,
      p,
      device = "pdf",
      dpi = 100,
      width = 10,
      height = 5
    )
  }
  p
}

#' @title Tables for KS test
#' @description Create tables for KS test
#' @param ks_results A data frame containing the results of the KS test when simulated from a single prior

KS_table <- function(ks_results) {
  # Remove duplicate rows based on the columns we want to keep
  ks_results_cleaned <- ks_results %>%
    distinct(parameter, prior, p_value, .keep_all = TRUE)

  # Escape underscores in the parameter column
  ks_results_cleaned$parameter <- gsub("_", "\\\\_", ks_results_cleaned$parameter)

  # Replace prior names with LaTeX-friendly names
  ks_results_cleaned$prior <- recode(ks_results_cleaned$prior,
    pc = "$\\pi_{\\mathrm{PC}}$",
    EG = "$\\pi_{\\mathrm{EG}}$",
    uniform = "$\\pi_{\\mathrm{U}}$",
    beta = "$\\pi_{\\beta}$"
  )

  # Format the p-values in scientific notation
  ks_results_cleaned$p_value <- format(ks_results_cleaned$p_value, scientific = TRUE, digits = 3)

  # Subset the data frame to only include the necessary columns for the table
  pvals_table <- ks_results_cleaned[, c("prior", "parameter", "p_value")] %>%
    pivot_wider(names_from = parameter, values_from = p_value)

  # Convert column names to LaTeX-friendly format
  colnames(pvals_table) <- c("Prior", "$\\log(\\kappa)$", "$v_1$", "$v_2$", "$\\log(\\sigma_u)$", "$\\log(\\sigma_{\\epsilon})$")

  # Generate LaTeX code for the table
  print(xtable(pvals_table, align = "lrrrrrr", digits = 3),
    sanitize.text.function = function(x) x,
    include.colnames = TRUE,
    include.rownames = FALSE
  )
}
