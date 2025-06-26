# Load the necessary packages
library(foreach)
library(doParallel)
library(microbenchmark)

# Register a parallel backend using all available cores
registerDoParallel(cores = detectCores())

# Define the number of iterations and the size of the random number vector to generate in each iteration
n_iterations <- 1000
size <- 10000

# Create a benchmark
benchmark_result <- microbenchmark(
  # Generate random numbers using a standard for loop
  "for loop" = {
    results <- numeric(n_iterations)
    for (i in 1:n_iterations) {
      results[i] <- sum(rnorm(size))
    }
  },

  # Generate random numbers using a parallel foreach loop
  "foreach" = {
    results <- foreach(i = 1:n_iterations, .combine = 'c') %dopar% {
      sum(rnorm(size))
    }
  },

  times = 10
)

# Print the benchmark result
print(benchmark_result)
