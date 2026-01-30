# Reproducibility Materials for Submission

This repository contains the R package `SPDEaniso` and all code necessary to reproduce the computational results for the manuscript "A Parameterization of Anisotropic Gaussian Fields with Penalized Complexity Priors".

The package is used to simulate from and perform Bayesian inference on anisotropic Gaussian fields using the anisotropic stochastic partial differential equations (SPDEs):

$$(\kappa^2-\nabla\cdot \mathbf{H}({\mathbf{v}})\nabla)u=\sqrt{4\pi}\kappa\sigma\mathcal{W}$$

where the parameters $\theta:=(\kappa, \mathbf{v}, \sigma)$ control the length scale, anisotropy, and variance of the field.

## System Requirements

- **R version:** R `(>= 3.6)` is required. R version `4.2.0` or newer is recommended.
- **Git:** Required for cloning the repository.
- **External Dependencies:** Some R packages used in the analysis may require external system libraries.
  - The `gsl` package requires the GNU Scientific Library. On Debian/Ubuntu, this can be installed with `sudo apt-get install libgsl-dev`.
  - Packages in the `INLA`/`fmesher` ecosystem may require a Fortran compiler (`gfortran`).

## Workflow for Reproducibility

Please follow these steps in order. The process involves first downloading all necessary files, then installing the required R packages, and finally running the analysis scripts.

### Step 1: Clone the Repository

First, clone this repository to your local machine using `git`. This will create a folder named `Aniso` containing all the required R scripts and package source code.

```bash
git clone https://github.com/AnonymousScholar2025/Aniso.git
```

### Step 2: Prepare the R Session

All subsequent commands must be run from within an R session that has the root of the cloned `Aniso` directory as its working directory.

- **Recommended Method:** This repository contains an RStudio Project file (`Aniso.Rproj`). Simply open this file in RStudio, which will automatically start an R session in the correct directory.
- **Alternative Method:** If not using RStudio, start a new R session and manually set the working directory before proceeding. For example:
  ```r
  # The path will depend on where you cloned the repository
  setwd("path/to/your/Aniso")
  ```

### Step 3: Install All R Packages

Once your R session is running in the correct directory, run the following commands to install all necessary packages. The full installation may take 5-20 minutes.

```r
# Step A: Install prerequisite and analysis packages from CRAN
# These packages are required to install the main package or to run the analysis scripts.
install.packages(c(
  "remotes",
  "devtools",
  "latex2exp",
  "lamW",
  "patchwork",
  "patchwork",
  "reshape2",
  "pracma",
  "gridExtra",
  "cowplot"
))

# Step B: Install the INLA ecosystem packages if they are missing
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(c("INLA", "inlabru", "fmesher"),
                   repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"),
                   dep = TRUE)
}

# Step C: Install the 'SPDEaniso' package from this repository
# This will automatically find and install your custom fmesher from the Remotes field
# in the DESCRIPTION file.
remotes::install_github("AnonymousScholar2025/Aniso", subdir = "Package")
# OR, if running locally:
# devtools::install("Package") # Must be run from the root directory of the repo
```

### Step 4: Run Analysis Scripts

After all packages are successfully installed, you can reproduce the results from the manuscript by running the scripts located in the `Manuscript/` directory.

> **Note:** If you want to use `SPDEaniso` functions interactively (outside of the provided scripts), you must first load the package with `library(SPDEaniso)`.

**Quick Examples (run in seconds):**

```r
# Section 3: Covariance plots
source("Manuscript/Section_3/Stationary_u_covariance_plot.R")
source("Manuscript/Section_3/Nonstationary_u_samples_and_covariance.R")

# Section 4: PC prior visualization
source("Manuscript/Section_4/PC_prior_plots_kappa_v.R")

# Section 5: Prior/posterior comparison plots (uses pre-computed data)
source("Manuscript/Section_5/Simulation_plots/Prior_Posterior_plots.R")
source("Manuscript/Section_5/Simulation_plots/Plot_map_loop_results.R")

# Section 6: Scores and prior comparison plots (uses pre-computed data)
source("Manuscript/Section_6/Plotting/Scores_sim_plot.R")
source("Manuscript/Section_6/Plotting/Prior_comparison_precip_rho_r.R")
source("Manuscript/Section_6/Plotting/Scores_Plots.R")  # Score comparison across observation counts
```

**Long-Running Simulations (hours):**

```r
# Section 5: Main simulation loop (600 iterations by default, takes hours)
source("Manuscript/Section_5/map_loop.R")

# Section 6: Precipitation simulations (100 loops by default, takes hours)
source("Manuscript/Section_6/Simulation_precip.R")
source("Manuscript/Section_6/Precipitation_loop.R")
```

### Directory Structure

- `Manuscript/Section_5/`: Contains code for defining and testing the PC priors.
  - `Prior_Posterior_Data/`: Contains calculated data for comparing the various priors and posteriors.
    > **Note:** The data in this folder is loaded by default in `Simulation_plots/Prior_Posterior_plots.R` to save time. You can regenerate it by running the `prior_posterior_plotter` functions within that script (warning: this may take time).
  - `map_loop.R`: The main simulation script for Section 5 results.
    > **Warning:** This script runs 600 iterations by default and is computationally intensive (hours). To test it quickly, we recommend reducing `number_of_loops` (line 116) to a small number (e.g., 2).
  - `Simulation_Results/`: Where the output of `map_loop.R` is saved.
  - `Simulation_plots/Plot_map_loop_results.R`: Generates the plots seen in the manuscript using the data from `Simulation_Results`.

- `Manuscript/Section_6/`: Contains code for the precipitation application (Section 6).
  - `Simulation_precip.R`: The main script for generating the simulation data.
    > **Note:** Data is saved to `Precip/data_sim/` and `Precip/Matlab/sim_data/`.
  - `Precip/data_sim/`: Stores the simulated precipitation data.
  - `Plotting/`: Contains scripts for reproducing manuscript plots.
    - `Scores_sim_plot.R`: Generates score comparison plots using the pre-calculated LOO results in `Precip/Results/`.
- `Manuscript/Section_3/`: Contains code for the theoretical properties of the field.

### Viewing .nb Files

Files with the `.nb` extension are Mathematica Notebook files. They can be viewed using Wolfram Mathematica or the free Wolfram Player.
