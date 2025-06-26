# SPDEaniso: Tool for working with anisotropic fields using SPDEs

This repository contains the R package and code used for the manuscript "A Parameterization of Anisotropic Gaussian Fields with Penalized Complexity Priors". The package uses the `fmesher` package and can be used to simulate from and perform Bayesian inference on anisotropic Gaussian fields using the anisotropic stochastic partial differential equations (SPDEs):

```math
(\kappa^2-\nabla\cdot \mathbf{H}({\mathbf{v}})\nabla)u=\kappa\sigma\mathcal{W}.
```

The parameters are $\kappa >0$, a two-dimensional vector $\mathbf{v}$, and $\sigma>0$. These parameters control the length scale, anisotropy, and variance, respectively.

For Bayesian inference, we consider a linear, noisy observation process:

```math
\mathbf{y} = \mathbf{A}\mathbf{u} + \mathbf{\epsilon},
```

where $\mathbf{A}$ is an observation matrix and $\mathbf{\epsilon}\sim\mathcal{N}(0,\sigma_{\mathbf{\epsilon}}^2\mathbf{I})$ is a vector of independent Gaussian noise. The package supports Bayesian inference on the parameters $\theta:=(\kappa, \mathbf{v}, \sigma,\sigma_{\mathbf{\epsilon}})$.

## Repository Structure

*   **`R/`**: Contains the source code for all functions in the `SPDEaniso` package.
*   **`Manuscript/`**: Contains the R scripts used to generate the figures and results presented in the manuscript.
*   **`tests/`**: Contains unit tests for the package.
*   **`DESCRIPTION`**: Lists all package dependencies and metadata.

## Installation

To run the code and reproduce the analyses, you will need to install this R package and its dependencies.

1.  **Install `devtools`**:
    If you don't already have it, install `devtools` from CRAN:
    ````r
    install.packages("devtools")
    ````

2.  **Install `SPDEaniso` and dependencies**:
    Clone this repository to your local machine. Then, from the R console, with the working directory set to the root of this project, run the following commands. The first command installs all the required dependencies listed in the [`DESCRIPTION`](DESCRIPTION) file. The second command installs the `SPDEaniso` package itself.
    ````r
    # Install all dependencies
    devtools::install_deps(dependencies = TRUE)

    # Install the package
    devtools::install()
    ````

## Reproducing Manuscript Results

After successfully installing the package and its dependencies, you can reproduce the results from the manuscript by running the scripts located in the `Manuscript/` directory.

**Note:** Some analysis scripts, like those for the precipitation data, may automatically download required datasets if they are not found locally.

For example, to run the main simulation loop from Section 3, execute the following command in your R console:

````r
source("Manuscript/Section_3/map_loop.R")