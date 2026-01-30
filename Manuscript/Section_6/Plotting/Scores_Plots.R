# =============================================================================
# Scores_Plots.R
# =============================================================================
# DESCRIPTION: Generates leave-one-out score comparison plots.
#
# DATA REQUIREMENTS:
# This script requires pre-generated simulation data. Before running this script,
# you must first run:
#
#   source("Manuscript/Section_6/Precipitation_loop.R")
#
# This will generate the required .RData files in:
#   - Manuscript/Section_6/Precip/Results/
#
# WARNING: The Precipitation_loop.R script is computationally intensive and
# takes several hours to complete.
# =============================================================================
library(ggplot2)
n_locations <- c(233, 200, 150, 100, 50, 20)
n_weights <- 1000
n_loops <- 10
models <- c("pc", "EG", "iso")

# =============================================================================
# DATA AVAILABILITY CHECK
# =============================================================================
# Check if required data files exist (sample one file)
sample_data_path <- paste0(
  "Manuscript/Section_6/Precip/Results/LOO_pc_w=", n_weights,
  "_loc=", n_locations[1], "_n_loops=", n_loops, ".RData"
)

if (!file.exists(sample_data_path)) {
  stop(
    "Required simulation data not found!\n",
    "  Missing: ", sample_data_path, "\n\n",
    "  To generate the data, you must run Precipitation_loop.R MULTIPLE TIMES:\n\n",
    "  1. Open Manuscript/Section_6/Precipitation_loop.R\n",
    "  2. For EACH value of n_locations, set the variable (around line 69) and run:\n\n",
    "     n_locations: ", paste(n_locations, collapse = ", "), "\n\n",
    "     That is ", length(n_locations), " runs total (each takes several hours).\n\n",
    "  3. Example for one run:\n",
    "       n_locations <- ", n_locations[1], "\n",
    "       source('Manuscript/Section_6/Precipitation_loop.R')\n\n",
    "  WARNING: Each run is computationally intensive and takes several hours.\n"
  )
}
# =============================================================================

# Initialize lists to store the data
loo_error <- list()
CRPS <- list()
DSS <- list()

# Load the data for each model and each number of locations
for (model in models) {
  for (n in n_locations) {
    data <- readRDS(paste0("Manuscript/Section_6/Precip/Results/LOO_", model, "_w=", n_weights, "_loc=", n, "_n_loops=", n_loops, ".RData"))
    loo_error[[paste0(model, "_", n)]] <- data$loo_error
    CRPS[[paste0(model, "_", n)]] <- data$CRPS
    DSS[[paste0(model, "_", n)]] <- data$DSS
  }
}

# Create a data frame for each score
df_loo_error <- data.frame(
  model = rep(models, each = length(n_locations)),
  n_locations = rep(n_locations, times = length(models)),
  score = unlist(loo_error),
  score_type = "RMSE"
)

df_CRPS <- data.frame(
  model = rep(models, each = length(n_locations)),
  n_locations = rep(n_locations, times = length(models)),
  score = unlist(CRPS),
  score_type = "CRPS"
)

df_DSS <- data.frame(
  model = rep(models, each = length(n_locations)),
  n_locations = rep(n_locations, times = length(models)),
  score = unlist(DSS),
  score_type = "DSS"
)

# Combine the data frames
df <- rbind(df_loo_error, df_CRPS, df_DSS)


# Plot the scores
p <- ggplot(df, aes(x = factor(n_locations), y = score, color = model, group = model)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  labs(x = "Number of Observations", y = "Score") +
  facet_wrap(~score_type, scales = "free_y") +
  theme_minimal() +
  theme(text = element_text(size = 20))
print(p)
ggsave("Manuscript/Simulation_images/Precipitation/LO_scores.pdf", plot = p, width = 15, height = 5, dpi = 300)
