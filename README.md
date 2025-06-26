# Reproducibility Materials for JASA Manuscript

This repository contains the R package `SPDEaniso` and all code necessary to reproduce the computational results for the manuscript "A Parameterization of Anisotropic Gaussian Fields with Penalized Complexity Priors".

The package is used to simulate from and perform Bayesian inference on anisotropic Gaussian fields using the anisotropic stochastic partial differential equations (SPDEs):

$$(\kappa^2-\nabla\cdot \mathbf{H}({\mathbf{v}})\nabla)u=\kappa\sigma\mathcal{W}$$

where the parameters $\theta:=(\kappa, \mathbf{v}, \sigma)$ control the length scale, anisotropy, and variance of the field.

## System Requirements

* **R version:** R `(>= 3.6)` is required. R version `4.2.0` or newer is recommended.
* **Git:** Required for cloning the repository.
* **External Dependencies:** Some R packages used in the analysis may require external system libraries.
    * The `gsl` package requires the GNU Scientific Library. On Debian/Ubuntu, this can be installed with `sudo apt-get install libgsl-dev`.
    * Packages in the `INLA`/`fmesher` ecosystem may require a Fortran compiler (`gfortran`).

## Workflow for Reproducibility

Please follow these steps in order. The process involves first downloading all necessary files, then installing the custom R package, and finally running the analysis scripts.

### Step 1: Clone the Repository

First, clone this repository to your local machine using `git`. This will create a folder named `Aniso` containing all the required R scripts and package source code.

```bash
git clone [https://github.com/AnonymousScholar2025/Aniso.git](https://github.com/AnonymousScholar2025/Aniso.git)
```

### Step 2: Install the R Package and Dependencies

All subsequent commands should be run within an R session started from the root of the cloned `Aniso` directory.

Once your R session is running in the correct directory, run the following commands in order. The full installation may take 5-20 minutes, as `devtools` and `INLA` are large packages.

```r
# Step A: Install prerequisite packages
# The installation process requires 'devtools', and 'remotes' is used to install from GitHub.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Step B: Install the INLA package
# INLA is not on CRAN and must be installed from its own repository.
install.packages("INLA", repos = c(getOption("repos"), INLA = "[https://inla.r-inla-download.org/R/stable](https://inla.r-inla-download.org/R/stable)"), dep = TRUE)

# Step C: Install the 'SPDEaniso' package
# We install from the 'package/' subdirectory in this repository.
remotes::install_github("AnonymousScholar2025/Aniso", subdir = "package", force = TRUE)
```

### Step 3: Run Analysis Scripts

After the package is successfully installed, you can reproduce the results from the manuscript by running the scripts located in the `Manuscript/` directory.

**Note:** Some analysis scripts, like those for the precipitation data, may automatically download required datasets if they are not found locally.

For example, to run the main simulation loop from Section 3, execute the following command in your R console:

````r
source("Manuscript/Section_3/map_loop.R")