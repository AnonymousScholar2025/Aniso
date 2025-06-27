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

Please follow these steps in order. The process involves first downloading all necessary files, then installing the required R packages, and finally running the analysis scripts.

### Step 1: Clone the Repository

First, clone this repository to your local machine using `git`. This will create a folder named `Aniso` containing all the required R scripts and package source code.

```bash
git clone [https://github.com/AnonymousScholar2025/Aniso.git](https://github.com/AnonymousScholar2025/Aniso.git)
```

### Step 2: Install R Packages

#### Prerequisite: Set the Working Directory

All subsequent commands must be run within an R session that has the root of the cloned `Aniso` directory as its working directory.

* **Recommended Method:** This repository should contain an RStudio Project file (`Aniso.Rproj`). Simply open this file in RStudio, which will automatically start an R session in the correct directory.
* **Alternative Method:** If not using RStudio, start an R session and manually set the working directory. For example:
    ```r
    # The path will depend on where you cloned the repository
    setwd("path/to/your/Aniso")
    ```

#### Installation Commands

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
  "reshape2"
  # Add any other CRAN packages your analysis scripts require here
))

# Step B: Install the INLA ecosystem packages if they are missing
if (!requireNamespace("INLA", quietly = TRUE)) {
  install.packages(c("INLA", "inlabru", "fmesher"), 
                   repos = c(getOption("repos"), INLA = "[https://inla.r-inla-download.org/R/stable](https://inla.r-inla-download.org/R/stable)"), 
                   dep = TRUE)
}

# Step C: Install the 'SPDEaniso' package from this repository
# This will automatically find and install your custom fmesher from the Remotes field
# in the DESCRIPTION file.
remotes::install_github("AnonymousScholar2025/Aniso", subdir = "package")
```

### Step 3: Run Analysis Scripts

After all packages are successfully installed, you can reproduce the results from the manuscript by running the scripts located in the `Manuscript/` directory.

For example, this script generates the main simulation maps for Section 3 of the manuscript.

```r
source("Manuscript/Section_3/map_loop.R")
```