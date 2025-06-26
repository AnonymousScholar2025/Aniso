library(ggplot2)
library(patchwork)
library(devtools)
document()

# Inspect parameters used for data generation ---------------
LOO_pc <- readRDS("Manuscript/Section_6/Precip/Results/LOO_pc_20000.RData")
LOO_iso <- readRDS("Manuscript/Section_6/Precip/Results/LOO_iso_20000.RData")
true_params <- list(
  ANISO = LOO_pc$MAP$par,
  ISO = c(LOO_iso$MAP$par[1], 0, 0, LOO_iso$MAP$par[2], LOO_iso$MAP$par[3])
)
# LOO_sim_pc_big <- readRDS("Manuscript/Section_6/Precip/Results/LOO_pc_20000.RData")
# convert_to_range_and_non_log(t(LOO_sim_pc_big$MAP$par))
# LOO_sim_pc_big$m_ub__y[1:2]
# LOO_sim_iso_big <- readRDS("Manuscript/Section_6/Precip/data/LOO_iso_20000.RData")
# convert_to_range_and_non_log(t(LOO_sim_iso_big$MAP$par))
# LOO_sim_iso_big$m_ub__y[1:2]

# Load the data for each Model and each number of locations ---------------
n_loc_ANISO <- c(2000, 1000, 800, 600, 400, 200, 100, 50, 25)
n_loc_ISO <- c(2000, 1000, 800, 600, 400, 200, 100, 50, 25)
# n_loc_ANISO <- c(50, 25)
# n_loc_ISO <- c(50, 25)
n_weights <- 100
n_loops <- 100
loo_error <- list()
CRPS <- list()
DSS <- list()
CIScore <- list()
map_dist_list <- list()
sim_types <- c("ANISO", "ISO")
Models <- c("pc", "EG", "iso")

# Load the data for each Model and each number of locations
for (sim_type in sim_types) {
  n_loc <- if (sim_type == "ANISO") {
    n_loc_ANISO
  } else {
    n_loc_ISO
  }
  if (sim_type == "ANISO") {
    sim_type_string <- "ANISO"
  } else if (sim_type == "ISO") {
    sim_type_string <- "ISO"
  }
  for (Model in Models) {
    for (n in n_loc) {
      data <- readRDS(
        paste0(
          "Manuscript/Section_6/Precip/data_sim/",
          sim_type_string,
          "/LOO_",
          Model,
          "_w=",
          n_weights,
          "_loc=",
          n,
          "_n_loops=",
          n_loops,
          ".RData"
        )
      )
      loo_error[[paste0(Model, "_", sim_type, "_", n)]] <- data$loo_error
      CRPS[[paste0(Model, "_", sim_type, "_", n)]] <- data$CRPS
      DSS[[paste0(Model, "_", sim_type, "_", n)]] <- data$DSS
      CI_value <- data$interval_scores_parameters
      if (Model == "iso") {
        CI_value <- c(CI_value[1], 0, 0, CI_value[2], CI_value[3])
      }
      names(CI_value) <- c("log_kappa", "v1", "v2", "log_sigma_u", "log_sigma_epsilon")
      CIScore[[paste0(Model, "_", sim_type, "_", n)]] <- CI_value
      if (Model == "pc" || Model == "EG") {
        diff <- data$MAP$par - true_params[[sim_type]]
        map_dist_list[[paste0(Model, "_", sim_type, "_", n)]] <- sqrt(sum(diff^
          2))
      } else if (Model == "iso") {
        diff <- c(data$MAP$par[1], 0, 0, data$MAP$par[2], data$MAP$par[3]) - true_params[[sim_type]]
        map_dist_list[[paste0(Model, "_", sim_type, "_", n)]] <- sqrt(sum(diff^
          2))
      }
    }
  }
}


# Data frame for each score -----------------------------------------------
df_loo_error <- data.frame(
  sim_type = rep(sim_types, each = length(Models) * length(n_loc)),
  Model = rep(Models, each = length(n_loc)),
  n_loc = rep(n_loc, times = length(Models) * length(sim_types)),
  score = unlist(loo_error),
  score_type = "RMSE"
)

df_CRPS <- data.frame(
  sim_type = rep(sim_types, each = length(Models) * length(n_loc)),
  Model = rep(Models, each = length(n_loc)),
  n_loc = rep(n_loc, times = length(Models) * length(sim_types)),
  score = unlist(CRPS),
  score_type = "CRPS"
)

df_DSS <- data.frame(
  sim_type = rep(sim_types, each = length(Models) * length(n_loc)),
  Model = rep(Models, each = length(n_loc)),
  n_loc = rep(n_loc, times = length(Models) * length(sim_types)),
  score = unlist(DSS),
  score_type = "DSS"
)
# Convert CIScore list to data frame properly
df_CI <- do.call(rbind, lapply(names(CIScore), function(name) {
  sim_type <- ifelse(grepl("ANISO", name), "ANISO", "ISO")
  model <- sub("_.*", "", name)
  n_loc <- as.numeric(sub(".*_", "", name))
  parameters <- CIScore[[name]]
  # Ensure parameters are named correctly
  data.frame(
    sim_type = rep(sim_type, length(parameters)),
    Model = rep(model, length(parameters)),
    n_loc = rep(n_loc, length(parameters)),
    parameter = names(parameters),
    value = unlist(parameters),
    stringsAsFactors = FALSE
  )
}))
rownames(df_CI) <- NULL
df_CI$parameter <- factor(df_CI$parameter, levels = c("log_kappa", "v1", "v2", "log_sigma_u", "log_sigma_epsilon"))
df_CI[df_CI == 0] <- NA
df_CI_filtered <- subset(df_CI, parameter != "log_sigma_epsilon")
latex_labels <- c(
  log_kappa = "log(kappa)",
  v1 = "v[1]",
  v2 = "v[2]",
  log_sigma_u = "log(sigma[u])",
  log_sigma_epsilon = "log(sigma[epsilon])"
)
txt_size <- 12
p_CI <- ggplot(df_CI, aes(x = n_loc, y = value, color = Model, group = Model)) +
  geom_line() +
  geom_point() +
  facet_grid(sim_type ~ parameter, scales = "free_y",
            #labeller = labeller(parameter = as_labeller(latex_labels, label_parsed))
            ) +
  labs(
    x = "Number of Locations",
    y = "CI score",
    color = "Model"
  ) +
  theme_minimal() +
  theme(text = element_text(size = 1.25 * txt_size))
p_CI
ggsave(
  "Manuscript/Simulation_images/Precipitation/Simulation_precip_CI_score.pdf",
  p_CI,
  width = 10,
  height = 6,
  dpi = 300
)



# Combine the data frames
df <- rbind(df_loo_error, df_CRPS, df_DSS)

#
p <- ggplot(df, aes(
  x = factor(n_loc),
  y = score,
  color = Model,
  group = Model
)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Observations", y = "Score") +
  facet_grid(sim_type ~ score_type, scales = "free_y") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12))
p

# Split the data frame by sim_type and score_type
df_split <- split(df, list(df$sim_type, df$score_type))

# Plot of scores ----------------------------------------------------------
for (name in names(df_split)) {
  # Extract the sim_type and score_type from the name
  sim_type <- strsplit(name, "\\.")[[1]][1]
  score_type <- strsplit(name, "\\.")[[1]][2]

  # Determine the full title based on sim_type
  sim_type_full <- ifelse(sim_type == "ANISO", "Anisotropic Data", "Isotropic Data")
  full_title <- paste0(sim_type_full, ": ", score_type)

  # Create the plot
  p <- ggplot(df_split[[name]], aes(
    x = factor(n_loc),
    y = score,
    color = Model,
    group = Model
  )) +
    geom_line() +
    geom_point() +
    labs(x = "Number of Observations", y = "Score", title = full_title) +
    theme_minimal() +
    theme(
      text = element_text(size = txt_size),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  # Save the plot as p1, p2, ..., p6
  assign(paste0("p", which(names(df_split) == name)), p, envir = .GlobalEnv)
}

p_combined <- p1 + p3 + p5 + p2 + p4 + p6 + plot_layout(guides = "collect")
print(p_combined)
ggsave(
  "Manuscript/Simulation_images/Precipitation/Simulation_precip_scores_combined.pdf",
  p_combined,
  width = 10,
  height = 6,
  dpi = 300
)



# Plot of difference of scores --------------------------------------------
df_diff <- df
pc_scores <- tapply(
  df_diff$score[df_diff$Model == "pc"],
  list(df_diff$sim_type[df_diff$Model == "pc"], df_diff$n_loc[df_diff$Model == "pc"], df_diff$score_type[df_diff$Model == "pc"]),
  mean
)
EG_scores <- tapply(
  df_diff$score[df_diff$Model == "EG"],
  list(df_diff$sim_type[df_diff$Model == "EG"], df_diff$n_loc[df_diff$Model == "EG"], df_diff$score_type[df_diff$Model == "EG"]),
  mean
)
iso_scores <- tapply(
  df_diff$score[df_diff$Model == "iso"],
  list(df_diff$sim_type[df_diff$Model == "iso"], df_diff$n_loc[df_diff$Model == "iso"], df_diff$score_type[df_diff$Model == "iso"]),
  mean
)


# Calculate the difference of each score with respect to the pc score
EG_diff_scores <- EG_scores - pc_scores
iso_diff_scores <- iso_scores - pc_scores

# Convert the arrays to data frames
EG_diff_df <- as.data.frame.table(EG_diff_scores)
iso_diff_df <- as.data.frame.table(iso_diff_scores)

# Rename the columns
names(EG_diff_df) <- c("sim_type", "n_loc", "score_type", "score")
names(iso_diff_df) <- c("sim_type", "n_loc", "score_type", "score")

# Add a Model column
EG_diff_df$Model <- "EG"
iso_diff_df$Model <- "iso"

# Combine the data frames
df_diff <- rbind(iso_diff_df, EG_diff_df)

# Create a plot for each subset
df_diff_split <- split(df_diff, list(df_diff$sim_type, df_diff$score_type))
for (name in names(df_diff_split)) {
  # Extract the sim_type and score_type from the name
  sim_type <- strsplit(name, "\\.")[[1]][1]
  score_type <- strsplit(name, "\\.")[[1]][2]

  # Create the plot
  d <- ggplot(
    df_diff_split[[name]],
    aes(
      x = factor(n_loc),
      y = score,
      color = Model,
      group = Model
    )
  ) +
    geom_line() +
    geom_point() +
    labs(
      x = "Number of Observations",
      y = "Score Difference",
      title = paste(sim_type, score_type)
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = txt_size),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_color_discrete(labels = c("EG-PC", "ISO-PC"))

  # Save the plot as d1, d2, ..., d6
  assign(paste0("d", which(names(df_diff_split) == name)), d, envir = .GlobalEnv)
}

d_combined <- d1 + d3 + d5 + d2 + d4 + d6 + plot_layout(guides = "collect")
print(d_combined)
custom_labeller <- function(variable, value) {
  if (variable == "sim_type") {
    ifelse(as.character(value) == "ANISO",
      "Anisotropic Data",
      "Isotropic Data"
    )
  } else {
    return(as.character(value))
  }
}
d <- ggplot(df_diff, aes(
  x = factor(n_loc),
  y = score,
  color = Model,
  group = Model
)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Observations", y = "Score Difference") +
  facet_grid(sim_type ~ score_type, scales = "free", labeller = custom_labeller) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = txt_size),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_discrete(labels = c("EG-PC", "ISO-PC"))

print(d)
ggsave(
  "Manuscript/Simulation_images/Precipitation/Simulation_precip_scores_diff_combined.pdf",
  d,
  width = 10,
  height = 6,
  dpi = 300
)


# v_abs <- list()
# for (sim_type in sim_types) {
#   n_loc <- if (sim_type == "ANISO") {
#     n_loc_ANISO
#   } else {
#     n_loc_ISO
#   }
#   for (Model in Models[1:2]) {
#     for (n in n_loc) {
#       data <- readRDS(
#         paste0(
#           "Precip/data_sim/",
#           sim_type_string,
#           "/LOO_",
#           Model,
#           "_w=",
#           n_weights,
#           "_loc=",
#           n,
#           "_n_loops=",
#           n_loops,
#           ".RData"
#         )
#       )
#       v_abs[[paste0(Model, "_", sim_type, "_", n)]] <- sqrt(sum(data$MAP_non_log_par[2:3] ^
#                                                                   2))
#     }
#   }
# }

# df_v_abs <- data.frame(
#   sim_type = rep(sim_types, each = length(Models[1:2]) * length(n_loc)),
#   Model = rep(Models[1:2], each = length(n_loc)),
#   n_loc = rep(n_loc, times = length(Models[1:2]) * length(sim_types)),
#   v_abs = unlist(v_abs)
# )
# df_v_abs_split <- split(df_v_abs, list(df_v_abs$sim_type, df_v_abs$Model))
#
# plots <- list()
#
# # Create a plot for each sim_type
# for (sim_type in sim_types) {
#   # Filter the data for the current sim_type
#   df_sim_type <- df_v_abs[df_v_abs$sim_type == sim_type, ]
#
#   # Create the plot
#   p <- ggplot(df_sim_type,
#               aes(
#                 x = factor(n_loc),
#                 y = v_abs,
#                 color = Model,
#                 group = Model
#               )) +
#     geom_line() +
#     geom_point() +
#     labs(x = "Number of Observations", y = "v_abs", title = sim_type) +
#     theme_minimal()
#
#   # Add the plot to the list
#   plots[[sim_type]] <- p
# }
#
# # Combine the plots
# p_combined <- plots[["ISO"]] + plots[["ANISO"]] +
#   plot_layout(guides = "collect") +
#   theme(strip.text = element_text(size = 12))
#
#
# # Display the combined plot
# print(p_combined)
#
# # Save the combined plot
# ggsave(
#   "Manuscript/Simulation_images/Simulation_precip_v_abs_combined.pdf",
#   p_combined,
#   width = 10,
#   height = 6,
#   dpi = 300
# )
#
