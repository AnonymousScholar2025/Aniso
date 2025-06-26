library(ggplot2)
n_locations <- c(233, 200, 150, 100, 50, 20)
n_weights <- 1000
n_loops <- 10
models <- c("pc", "EG", "iso")

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
ggplot(df, aes(x = factor(n_locations), y = score, color = model, group = model)) +
  geom_line() +
  geom_point() +
  labs(x = "Number of Observations", y = "Score") +
  facet_wrap(~score_type, scales = "free_y") +
  theme_minimal() +
  theme(text = element_text(size = 20))
ggsave("Manuscript/Simulation_images/Precipitation/LO_scores.pdf", width = 15, height = 5, dpi = 300)
