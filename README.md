# Reproducibility Materials for JASA Manuscript

This repository contains the R package `SPDEaniso` and all code necessary to reproduce the computational results for the manuscript "A Parameterization of Anisotropic Gaussian Fields with Penalized Complexity Priors".

The package is used to simulate from and perform Bayesian inference on anisotropic Gaussian fields using the anisotropic stochastic partial differential equations (SPDEs):

$$(\kappa^2-\nabla\cdot \mathbf{H}({\mathbf{v}})\nabla)u=\kappa\sigma\mathcal{W}$$

where the parameters $\theta:=(\kappa, \mathbf{v}, \sigma)$ control the length scale, anisotropy, and variance of the field.

## System Requirements

* **R version:** R `(>= 3.6)` is required. R version `4.2.0` or newer is recommended.
* **External Dependencies:** Some R packages used in the analysis may require external system libraries.
    * The `gsl` package requires the GNU Scientific Library. On Debian/Ubuntu, this can be installed with `sudo apt-get install libgsl-dev`.
    * Packages in the `INLA`/`fmesher` ecosystem may require a Fortran compiler (`gfortran`).

## Installation

The `SPDEaniso` R package and all of its dependencies can be installed from this GitHub repository using the `remotes` package. This command will automatically install all required packages from CRAN and the correct `fmesher` branch from GitHub as specified in the `DESCRIPTION` file.

```r
# If 'remotes' is not installed, run this first:
# install.packages("remotes")

# Install the package and all dependencies from this repository
remotes::install_github("AnonymousScholar2025/JASA_Aniso_Submission_2")
```

## Reproducing Manuscript Results

All scripts to reproduce the figures and tables in the manuscript are located in the `Manuscript/` directory. Please run scripts from the project's root directory.


### To Reproduce All Results
After successfully installing the package and its dependencies, you can reproduce the results from the manuscript by running the scripts located in the `Manuscript/` directory.

**Note:** Some analysis scripts, like those for the precipitation data, may automatically download required datasets if they are not found locally.

For example, to run the main simulation loop from Section 3, execute the following command in your R console:

````r
source("Manuscript/Section_3/map_loop.R")